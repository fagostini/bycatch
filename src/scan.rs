//! Full pass over the FASTQ, classify every Nth read (or read pair).
use crate::index::{Index, SHARED};
use crate::kmer::{kmers, low_complexity, minimizers, mix, K};
use crate::Res;
use needletail::{parse_fastx_file, parse_fastx_stdin, FastxReader};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::sync::OnceLock;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread::JoinHandle;

fn open(path: &str) -> Res<Box<dyn FastxReader>> {
    let r = if path == "-" { parse_fastx_stdin() } else { parse_fastx_file(path) };
    r.map_err(|e| format!("{path}: {e}").into())
}

/// Below this diversity an organism's hits pile on few markers: a conserved gene
/// (e.g. 16S rRNA from an organism not in the index), a repeat or an amplicon.
/// Real genomes measured 0.59 to 1.0; 16S artefacts 0.005 to 0.16.
const LOW_DIVERSITY: f64 = 0.25;
/// Too few hits to judge.
const MIN_HITS_FOR_FLAG: u64 = 20;

/// Marker density of an organism's assigned fragments against what reads copied from its
/// reference would show (its unique fraction), and the per-base identity that ratio
/// suggests: density ≈ identity^k when every k-mer needs k matching bases. Sequencing
/// errors lower it too, so it is a relative measure, not an ANI.
fn divergence(own_hits: u64, minimizers: u64, unique_fraction: f64) -> Option<(f64, f64)> {
    if own_hits < MIN_HITS_FOR_FLAG || minimizers == 0 || unique_fraction <= 0.0 {
        return None;
    }
    let density = own_hits as f64 / minimizers as f64 / unique_fraction;
    Some((density, density.min(1.0).powf(1.0 / K as f64)))
}

/// Distinct markers hit, divided by the number expected if the same number of
/// hits fell uniformly over the organism's markers: markers * (1 - e^(-hits/markers)).
fn diversity(hits: u64, distinct: u64, markers: u64) -> Option<f64> {
    if hits == 0 || markers == 0 {
        return None;
    }
    let m = markers as f64;
    Some(distinct as f64 / (-m * (-(hits as f64) / m).exp_m1()))
}

/// A sampled R2 record: its sequence, plus id and qualities when reads are written out.
type Mate = (Vec<u8>, Option<Raw>);

/// Reads R2 on a background thread and hands over the sequences of the sampled records.
struct MateReader {
    rx: Receiver<Vec<Mate>>,
    batch: std::vec::IntoIter<Mate>,
    thread: Option<JoinHandle<Res<u64>>>,
}

impl MateReader {
    fn start(path: String, every: u64, raw: bool) -> MateReader {
        let (tx, rx) = sync_channel::<Vec<Mate>>(16);
        let thread = std::thread::spawn(move || -> Res<u64> {
            let mut reader = open(&path)?;
            let (mut n, mut batch) = (0u64, Vec::new());
            while let Some(rec) = reader.next() {
                let rec = rec.map_err(|e| format!("{path}: {e}"))?;
                n += 1;
                if (n - 1).is_multiple_of(every) {
                    let r = if raw { Some(Raw::of(&rec, &path)?) } else { None };
                    batch.push((rec.seq().into_owned(), r));
                    if batch.len() == 1024 && tx.send(std::mem::take(&mut batch)).is_err() {
                        break; // main thread stopped on an error
                    }
                }
            }
            let _ = tx.send(batch);
            Ok(n)
        });
        MateReader { rx, batch: Vec::new().into_iter(), thread: Some(thread) }
    }

    fn next(&mut self) -> Res<Mate> {
        loop {
            if let Some(s) = self.batch.next() {
                return Ok(s);
            }
            match self.rx.recv() {
                Ok(b) => self.batch = b.into_iter(),
                // thread ended: report its error, else R2 ran out first
                Err(_) => {
                    let n2 = self.finish_thread()?;
                    return Err(format!("R2 has fewer reads than R1 (R2: {n2})").into());
                }
            }
        }
    }

    fn finish_thread(&mut self) -> Res<u64> {
        let t = self.thread.take().ok_or("R2 reader already finished")?;
        t.join().map_err(|_| "R2 reader thread panicked")?
    }

    /// Total R2 records, once R2 is fully read.
    fn finish(mut self) -> Res<u64> {
        while self.rx.recv().is_ok() {} // drain, so the thread can reach the end of R2
        self.finish_thread()
    }

    /// Stop reading R2 before its end; reports only a read error.
    fn stop(self) -> Res<()> {
        let MateReader { rx, thread, .. } = self;
        drop(rx); // the thread's next send fails and it returns
        match thread {
            Some(t) => t.join().map_err(|_| "R2 reader thread panicked")?.map(|_| ()),
            None => Ok(()),
        }
    }
}

/// Sampled fragments handed from the reader thread to the classifier threads at once.
const BATCH: usize = 4096;

/// Adapters read through after a short insert, and primer dimers: Illumina TruSeq
/// (read 1 and 2, and the P5 side), Nextera (read 1 and 2, mosaic end), the P5 and P7
/// flow-cell ends that follow the index, MGI (read 1 and 2).
/// shortcut: built-in list and exact 17-mers only; add a --adapters FASTA if other kits
/// show up as `unmatched`.
const ADAPTERS: [&[u8]; 10] = [
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    b"CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
    b"CTGTCTCTTATACACATCTGACGCTGCCGACGA",
    b"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    b"ATCTCGTATGCCGTCTTCTGCTTG",
    b"AATGATACGGCGACCACCGAGATCTACAC",
    b"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    b"AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
];
/// Short enough to catch a dimer that starts mid-adapter and goes on into the index
/// (as few as 18 adapter bases seen on a NovaSeq run); a random 300 bp pair matches
/// one of the ~600 adapter 17-mers with probability ~1e-5, and only reads with no
/// other hit are checked.
pub const ADAPTER_K: usize = 17;

pub fn adapter_kmers() -> &'static HashSet<u64> {
    static SET: OnceLock<HashSet<u64>> = OnceLock::new();
    SET.get_or_init(|| {
        let mut set = HashSet::new();
        for a in ADAPTERS {
            kmers(a, ADAPTER_K, |k| {
                set.insert(k);
            });
        }
        set
    })
}

/// Rows after the organisms, in report order.
pub const SPECIAL: [&str; 6] = ["shared", "weak", "conflict", "adapter", "low_complexity", "unmatched"];
const C_SHARED: usize = 0;
const C_WEAK: usize = 1;
const C_CONFLICT: usize = 2;
const C_ADAPTER: usize = 3;
const C_LOW: usize = 4;
const C_UNMATCHED: usize = 5;

/// Id and qualities of a record, kept only to write the reads out again.
struct Raw {
    id: Vec<u8>,
    qual: Vec<u8>,
}

impl Raw {
    fn of(rec: &needletail::parser::SequenceRecord, path: &str) -> Res<Raw> {
        let qual = rec.qual().ok_or_else(|| format!("{path}: --depleted-out needs FASTQ, not FASTA"))?;
        Ok(Raw { id: rec.id().to_vec(), qual: qual.to_vec() })
    }

    /// The record as FASTQ, with sequence `seq`.
    fn fastq(&self, seq: &[u8], out: &mut Vec<u8>) {
        for part in [b"@", &self.id[..], b"\n", seq, b"\n+\n", &self.qual[..], b"\n"] {
            out.extend_from_slice(part);
        }
    }
}

/// One sampled read, or read pair.
struct Fragment {
    rec: u64,    // record number in R1, from 1
    id: Vec<u8>, // empty unless read ids are written out
    seq: Vec<u8>,
    mate: Option<Vec<u8>>,
    raw: Option<(Raw, Option<Raw>)>, // R1 and R2 id and qualities, for --depleted-out
}

/// Reads R1 (and R2 through a MateReader) and sends the sampled fragments in batches.
/// Stops after `max_reads` fragments, or when the receiver is dropped.
/// Returns the records read and whether it stopped before the end of the file.
fn read_fragments(
    r1: &str,
    r2: Option<&str>,
    every: u64,
    max_reads: Option<u64>,
    keep_ids: bool,
    keep_raw: bool,
    tx: SyncSender<Vec<Fragment>>,
) -> Res<(u64, bool)> {
    let mut reader = open(r1)?;
    let mut mates = r2.map(|p| MateReader::start(p.to_string(), every, keep_raw));
    let (mut n, mut sampled, mut stopped) = (0u64, 0u64, false);
    let mut batch = Vec::with_capacity(BATCH);
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| format!("{r1}: {e}"))?;
        n += 1;
        if !(n - 1).is_multiple_of(every) {
            continue;
        }
        let (mate, mate_raw) = match &mut mates {
            Some(m) => {
                let (seq, raw) = m.next()?;
                (Some(seq), raw)
            }
            None => (None, None),
        };
        let id = if keep_ids { rec.id().to_vec() } else { Vec::new() };
        let raw = if keep_raw { Some((Raw::of(&rec, r1)?, mate_raw)) } else { None };
        batch.push(Fragment { rec: n, id, seq: rec.seq().into_owned(), mate, raw });
        sampled += 1;
        if max_reads == Some(sampled) {
            stopped = true;
            break;
        }
        if batch.len() == BATCH && tx.send(std::mem::take(&mut batch)).is_err() {
            stopped = true; // the classifier stopped: --precision reached, or its own error
            break;
        }
    }
    let _ = tx.send(batch);
    if let Some(m) = mates {
        if stopped {
            m.stop()?;
        } else {
            let n2 = m.finish()?;
            if n2 != n {
                return Err(format!("R1 has {n} reads but R2 has {n2}").into());
            }
        }
    }
    Ok((n, stopped))
}

/// The strict rule. The top organism is assigned when it has at least 2 marker hits,
/// at least `min_frac` of the fragment's minimizers, and with its group's markers at
/// least 90 % of all marker hits. The `min_frac` test scales with read length: a long
/// host read can pick up 2 stray markers of another organism and none of its own.
/// Otherwise a group is assigned when its markers and its members' pass the same tests.
/// Returns the organism or group id and the hits that support it.
fn call(hits: &[(usize, u32, u32)], genomes: usize, group_of: &[Option<usize>], minimizers: usize, min_frac: f64) -> Option<(usize, u32)> {
    let total: u32 = hits.iter().map(|h| h.1).sum();
    let ok = |own: u32, support: u32| own >= 2 && support * 10 >= total * 9 && own as f64 >= min_frac * minimizers as f64;
    let count = |g: usize| hits.iter().find(|h| h.0 == g).map_or(0, |h| h.1);
    if let Some(&(g, own, _)) = hits.iter().filter(|h| h.0 < genomes).max_by_key(|h| h.1) {
        let support = own + group_of[g].map_or(0, count);
        if ok(own, support) {
            return Some((g, support));
        }
    }
    let in_group = |id: usize, gr: usize| id == gr || id < genomes && group_of[id] == Some(gr);
    let mut best: Option<(usize, u32)> = None;
    for gr in hits.iter().filter_map(|h| if h.0 < genomes { group_of[h.0] } else { Some(h.0) }) {
        let sum = hits.iter().filter(|h| in_group(h.0, gr)).map(|h| h.1).sum();
        if best.is_none_or(|b| sum > b.1) {
            best = Some((gr, sum));
        }
    }
    best.filter(|&(_, sum)| ok(sum, sum))
}

/// Index hits of one fragment.
struct Call {
    hits: Vec<(usize, u32, u32)>, // (organism, marker hits, markers not hit before in this scan)
    shared: u32,
    minimizers: usize,
    adapter: bool, // only looked for when nothing else matched
    key: Option<u64>, // start of the fragment, to find PCR and optical duplicates
}

/// Whether the read's first 27 bases can tell molecules apart: not a microsatellite
/// such as (GT)n, which many different molecules start with.
fn informative_start(seq: &[u8]) -> bool {
    let mut ok = false;
    kmers(&seq[..seq.len().min(K)], K, |k| ok = !low_complexity(k));
    ok
}

/// Hash of the first 32 bases of the read and of its mate: duplicates share both.
/// None unless both starts are informative: a pair whose R2 is poly-G (no signal)
/// would otherwise be keyed by R1 alone, and adapter dimers all share their R1.
fn dup_key(f: &Fragment) -> Option<u64> {
    let mate = f.mate.as_deref().unwrap_or(&[]);
    if !informative_start(&f.seq) || f.mate.is_some() && !informative_start(mate) {
        return None;
    }
    let bytes = f.seq.iter().take(32).chain(b"|").chain(mate.iter().take(32));
    Some(mix(bytes.fold(0xcbf29ce484222325u64, |h, &b| (h ^ b as u64).wrapping_mul(0x100000001b3))))
}

/// Per-thread scratch space of `classify`: minimizers, and their index ranges.
type Scratch = (Vec<u64>, Vec<(usize, usize)>);

fn classify(idx: &Index, seen: &[AtomicU64], f: &Fragment, (mins, ranges): &mut Scratch) -> Call {
    mins.clear();
    minimizers(&f.seq, idx.w, mins);
    if let Some(m) = &f.mate {
        minimizers(m, idx.w, mins);
    }
    mins.sort_unstable();
    mins.dedup();
    // Lookups are random memory reads: issue all of one stage before using any, so
    // the cache misses of a fragment overlap.
    mins.iter().for_each(|&k| idx.prefetch_table(k));
    ranges.clear();
    ranges.extend(mins.iter().map(|&k| idx.range(k)));
    ranges.iter().for_each(|r| idx.prefetch_entries(r.0));
    let (mut hits, mut shared) = (Vec::<(usize, u32, u32)>::new(), 0);
    for (&kmer, &(lo, hi)) in mins.iter().zip(ranges.iter()) {
        match idx.find_in(kmer, lo, hi) {
            None => {}
            Some((_, SHARED)) => shared += 1,
            Some((pos, id)) => {
                let bit = 1u64 << (pos % 64);
                // fetch_or returns the old word: exactly one thread sees the bit unset
                let new = (seen[pos / 64].fetch_or(bit, Relaxed) & bit == 0) as u32;
                match hits.iter_mut().find(|h| h.0 == id as usize) {
                    Some(h) => (h.1, h.2) = (h.1 + 1, h.2 + new),
                    None => hits.push((id as usize, 1, new)),
                }
            }
        }
    }
    hits.sort_unstable();
    let mut adapter = false;
    if hits.is_empty() && shared == 0 {
        let set = adapter_kmers();
        for s in std::iter::once(&f.seq).chain(&f.mate) {
            kmers(s, ADAPTER_K, |k| adapter |= set.contains(&k));
        }
    }
    Call { hits, shared, minimizers: mins.len(), adapter, key: dup_key(f) }
}

/// 95 % Wilson score interval for a fraction p of n trials, in percent. Stays inside
/// 0–100 and is usable at p = 0 (upper bound about 3/n, the "rule of three").
/// n may be fractional: the effective number of fragments after duplicates.
fn wilson(p: f64, n: f64) -> (f64, f64) {
    if n <= 0.0 {
        return (0.0, 100.0);
    }
    let z2 = 1.96f64 * 1.96;
    let centre = (p + z2 / (2.0 * n)) / (1.0 + z2 / n);
    let half = z2.sqrt() / (1.0 + z2 / n) * (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt();
    (100.0 * (centre - half).max(0.0), 100.0 * (centre + half).min(1.0))
}

/// Duplicate clusters among the sampled fragments. Duplicates are not independent, so
/// the CI uses n / D fragments, where D = Σm² / Σm over cluster sizes m (the mean
/// cluster size seen by a fragment; 1 without duplicates). Keeps the clusters of a
/// fixed hash fraction only, halved whenever more than MAX keys are held, so memory
/// stays bounded; the clusters kept are whole, so D stays unbiased.
struct Dups {
    shift: u32, // keys with their low `shift` bits zero are kept
    clusters: HashMap<u64, u32>,
    sum_m: u64,
    sum_m2: u64,
}

impl Dups {
    const MAX: usize = 1 << 21;

    fn new() -> Dups {
        Dups { shift: 0, clusters: HashMap::new(), sum_m: 0, sum_m2: 0 }
    }

    fn add(&mut self, key: u64) {
        if key & ((1u64 << self.shift) - 1) != 0 {
            return;
        }
        let m = self.clusters.entry(key).or_insert(0);
        self.sum_m2 += 2 * *m as u64 + 1; // (m+1)² - m²
        self.sum_m += 1;
        *m += 1;
        if self.clusters.len() > Self::MAX {
            self.shift += 1;
            let mask = (1u64 << self.shift) - 1;
            self.clusters.retain(|k, _| k & mask == 0);
            self.sum_m = self.clusters.values().map(|&m| m as u64).sum();
            self.sum_m2 = self.clusters.values().map(|&m| (m as u64).pow(2)).sum();
        }
    }

    fn design_effect(&self) -> f64 {
        (self.sum_m2 as f64 / self.sum_m.max(1) as f64).max(1.0)
    }
}

/// Splits the fragments the strict rule leaves as `weak` or `conflict` between the
/// organisms they hit, in proportion to abundance (EM, as Bracken does for Kraken).
/// Model: each marker hit of a fragment from organism g is a marker of g with
/// probability 1 - eps, else of another organism. Only organisms with at least 2 hits
/// in the fragment are candidates (the strict rule's 2-hit minimum), so 1-hit `weak`
/// reads stay unassigned. `ambiguous` maps a candidate pattern [(organism, hits)] to its
/// fragment count.
/// Returns expected fragments per organism, the strictly assigned ones included.
fn em(assigned: &[u64], ambiguous: &HashMap<Vec<(usize, u32)>, u64>, eps: f64) -> Vec<f64> {
    let total = (assigned.iter().sum::<u64>() + ambiguous.values().sum::<u64>()) as f64;
    let log_odds = ((1.0 - eps) / eps).ln();
    let mut pi = vec![1.0 / assigned.len() as f64; assigned.len()];
    let mut est: Vec<f64> = assigned.iter().map(|&x| x as f64).collect();
    if ambiguous.is_empty() {
        return est;
    }
    for _ in 0..1000 {
        est = assigned.iter().map(|&x| x as f64).collect();
        for (pattern, &n) in ambiguous {
            let top = pattern.iter().map(|p| p.1).max().unwrap_or(0) as f64;
            // likelihood ratio (1-eps)^c * eps^(T-c) ∝ e^(log_odds * c), scaled by the top count
            let w: Vec<f64> = pattern.iter().map(|&(g, c)| pi[g] * (log_odds * (c as f64 - top)).exp()).collect();
            let sum: f64 = w.iter().sum();
            if sum > 0.0 {
                for (&(g, _), w) in pattern.iter().zip(&w) {
                    est[g] += n as f64 * w / sum;
                }
            }
        }
        let change = est.iter().zip(&pi).map(|(e, p)| (e / total - p).abs()).fold(0.0, f64::max);
        pi = est.iter().map(|e| e / total).collect();
        if change < 1e-12 {
            break;
        }
    }
    est
}

pub struct Opts {
    pub every: u64,
    pub max_reads: Option<u64>,
    pub threads: usize,
    /// stop once the contamination CI's half-width is at most this many percentage points
    pub precision: Option<f64>,
    pub reads_out: Option<String>,
    pub unmatched_out: Option<String>,
    pub unmatched_n: usize,
    /// share of the fragment's minimizers the top organism must hit (0 = off)
    pub min_frac: f64,
    /// write the fragments not assigned to the host (or its group) to <prefix>[_R1/_R2].fastq.gz
    pub depleted_out: Option<String>,
}

/// Fragments sampled before --precision may stop the scan.
const MIN_FOR_PRECISION: u64 = 10_000;

/// Counts of one scanned sample.
pub struct Report {
    pub host: Option<usize>,
    /// the host's group: "host or a close relative", not contamination
    pub host_group: Option<usize>,
    pub reads: Vec<u64>,    // fragments assigned per organism
    pub hits: Vec<u64>,     // marker k-mer hits per organism
    pub distinct: Vec<u64>, // different markers hit per organism
    pub own_hits: Vec<u64>, // marker hits of each organism in its own assigned fragments
    pub own_minimizers: Vec<u64>, // minimizers of those fragments
    pub special: [u64; 6],  // fragments per SPECIAL category
    pub sampled: u64,
    /// records in the file, or up to the last sampled one when stopped early
    pub total: u64,
    pub stop: Option<&'static str>, // "max_reads" or "precision"
    pub design: f64,                // duplicate design effect, >= 1
    pub eps: f64,                   // per-hit error rate used by EM
    pub est: Vec<f64>,              // EM fragments per organism
}

impl Report {
    pub fn pct(&self, x: f64) -> f64 {
        100.0 * x / self.sampled.max(1) as f64
    }

    /// 95 % CI of x fragments, in percent, over the effective number of fragments.
    pub fn ci(&self, x: u64) -> (f64, f64) {
        wilson(x as f64 / self.sampled.max(1) as f64, self.sampled as f64 / self.design)
    }

    /// Whether fragments of organism or group `g` count as host.
    pub fn is_host(&self, g: usize) -> bool {
        self.host == Some(g) || self.host_group == Some(g)
    }

    /// Non-host fragments assigned by the strict rule (needs a host).
    pub fn contamination(&self) -> Option<u64> {
        self.host.map(|_| (0..self.reads.len()).filter(|&g| !self.is_host(g)).map(|g| self.reads[g]).sum())
    }

    pub fn em_contamination(&self) -> Option<f64> {
        self.host.map(|_| (0..self.est.len()).filter(|&g| !self.is_host(g)).map(|g| self.est[g]).sum())
    }

    pub fn write_tsv(&self, idx: &Index, out: &mut impl Write) -> Res<()> {
        let ci = |x: u64| {
            let (lo, hi) = self.ci(x);
            format!("{lo:.4}\t{hi:.4}")
        };
        writeln!(
            out,
            "category\treads\tpercent\tmarker_hits\tdistinct_markers\tindex_markers\tdiversity\tnote\tci95_low\tci95_high\tem_percent\tmarker_density\tidentity_est"
        )?;
        for g in 0..idx.names.len() {
            let div = diversity(self.hits[g], self.distinct[g], idx.markers[g]);
            let mut note = Vec::new();
            if self.host == Some(g) {
                note.push("host");
            }
            if g >= idx.genomes {
                note.push(if self.host_group == Some(g) { "group,host_group" } else { "group" });
            }
            if self.hits[g] >= MIN_HITS_FOR_FLAG && div.is_some_and(|d| d < LOW_DIVERSITY) {
                note.push("low_diversity");
            }
            let dv = (g < idx.genomes)
                .then(|| divergence(self.own_hits[g], self.own_minimizers[g], idx.unique_fraction(g)))
                .flatten();
            writeln!(
                out,
                "{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}",
                idx.names[g],
                self.reads[g],
                self.pct(self.reads[g] as f64),
                self.hits[g],
                self.distinct[g],
                idx.markers[g],
                div.map_or(String::new(), |d| format!("{d:.3}")),
                note.join(","),
                ci(self.reads[g]),
                self.pct(self.est[g]),
                dv.map_or(String::new(), |d| format!("{:.3}", d.0)),
                dv.map_or(String::new(), |d| format!("{:.4}", d.1))
            )?;
        }
        for (name, &x) in SPECIAL.iter().zip(&self.special) {
            writeln!(out, "{name}\t{x}\t{:.4}\t\t\t\t\t\t{}\t\t\t", self.pct(x as f64), ci(x))?;
        }
        writeln!(out, "# total_fragments\t{}", self.total)?;
        if let Some(s) = self.stop {
            writeln!(out, "# stopped_at_{s}\t{}", self.sampled)?;
        }
        writeln!(out, "# sampled_fragments\t{}", self.sampled)?;
        writeln!(out, "# duplicate_design_effect\t{:.3}", self.design)?;
        // an organism with no fragment is below this, at 95 % (the CI's upper end for 0)
        writeln!(out, "# absent_below_percent\t{:.4}", self.ci(0).1)?;
        writeln!(out, "# em_marker_error\t{:.6}", self.eps)?;
        if let (Some(c), Some(e)) = (self.contamination(), self.em_contamination()) {
            writeln!(out, "# contamination_percent\t{:.4}", self.pct(c as f64))?;
            writeln!(out, "# contamination_ci95\t{}", ci(c))?;
            writeln!(out, "# em_contamination_percent\t{:.4}", self.pct(e))?;
        }
        Ok(())
    }
}

pub fn host_id(idx: &Index, host: Option<&str>) -> Res<Option<usize>> {
    match host {
        Some(h) => Ok(Some(
            idx.names[..idx.genomes]
                .iter()
                .position(|n| n == h)
                .ok_or_else(|| format!("host '{h}' is not an organism in the index: {:?}", &idx.names[..idx.genomes]))?,
        )),
        None => Ok(None),
    }
}

/// (read id, sequence, mate sequence)
type Kept = (Vec<u8>, Vec<u8>, Option<Vec<u8>>);

/// Reservoir sample of `unmatched` fragments (Algorithm R, fixed seed: reproducible).
struct Reservoir {
    n: usize,
    seen: u64,
    rng: u64,
    kept: Vec<Kept>,
}

impl Reservoir {
    fn offer(&mut self, f: &Fragment) {
        self.seen += 1;
        let slot = if self.kept.len() < self.n {
            self.kept.push(Default::default());
            self.kept.len() - 1
        } else {
            self.rng = mix(self.rng.wrapping_add(0x9e3779b97f4a7c15));
            match (self.rng % self.seen) as usize {
                j if j < self.n => j,
                _ => return,
            }
        };
        self.kept[slot] = (f.id.clone(), f.seq.clone(), f.mate.clone());
    }

    fn write(&self, path: &str) -> Res<()> {
        let file = std::fs::File::create(path).map_err(|e| format!("{path}: {e}"))?;
        let mut w = std::io::BufWriter::new(file);
        for (id, seq, mate) in &self.kept {
            let id = String::from_utf8_lossy(id.split(|&b| b == b' ').next().unwrap_or(b""));
            match mate {
                Some(m) => {
                    writeln!(w, ">{id}/1\n{}", String::from_utf8_lossy(seq))?;
                    writeln!(w, ">{id}/2\n{}", String::from_utf8_lossy(m))?;
                }
                None => writeln!(w, ">{id}\n{}", String::from_utf8_lossy(seq))?,
            }
        }
        w.flush()?;
        Ok(())
    }
}

/// gzip writer of one FASTQ on its own thread, so compression runs beside classification.
struct GzThread {
    tx: SyncSender<Vec<u8>>,
    thread: JoinHandle<Res<()>>,
}

impl GzThread {
    fn start(path: String) -> Res<GzThread> {
        let file = std::fs::File::create(&path).map_err(|e| format!("{path}: {e}"))?;
        let (tx, rx) = sync_channel::<Vec<u8>>(4);
        let thread = std::thread::spawn(move || -> Res<()> {
            // level 1: output size matters less than keeping up with classification
            let mut w = flate2::write::GzEncoder::new(std::io::BufWriter::new(file), flate2::Compression::fast());
            for chunk in rx {
                w.write_all(&chunk).map_err(|e| format!("{path}: {e}"))?;
            }
            w.finish()?.flush().map_err(|e| format!("{path}: {e}"))?;
            Ok(())
        });
        Ok(GzThread { tx, thread })
    }

    fn send(&self, chunk: Vec<u8>) -> Res<()> {
        // a failed send means the writer stopped on an error: finish() reports it
        let _ = self.tx.send(chunk);
        Ok(())
    }

    fn finish(self) -> Res<()> {
        drop(self.tx);
        self.thread.join().map_err(|_| "--depleted-out writer thread panicked")?
    }
}

/// `--depleted-out`: the FASTQ of each batch is built in order, then handed to the writers.
struct Depleted {
    r1: Vec<u8>,
    r2: Vec<u8>,
    w1: GzThread,
    w2: Option<GzThread>,
}

impl Depleted {
    fn start(prefix: &str, paired: bool) -> Res<Depleted> {
        let (w1, w2) = if paired {
            (GzThread::start(format!("{prefix}_R1.fastq.gz"))?, Some(GzThread::start(format!("{prefix}_R2.fastq.gz"))?))
        } else {
            (GzThread::start(format!("{prefix}.fastq.gz"))?, None)
        };
        Ok(Depleted { r1: Vec::new(), r2: Vec::new(), w1, w2 })
    }

    fn flush_batch(&mut self) -> Res<()> {
        self.w1.send(std::mem::take(&mut self.r1))?;
        if let Some(w2) = &self.w2 {
            w2.send(std::mem::take(&mut self.r2))?;
        }
        Ok(())
    }

    fn finish(mut self) -> Res<()> {
        self.flush_batch()?;
        self.w1.finish()?;
        self.w2.map_or(Ok(()), GzThread::finish)
    }
}

/// Classifies one sample (a file, or a read pair).
pub fn scan_sample(idx: &Index, opts: &Opts, host: Option<usize>, r1: &str, r2: Option<&str>) -> Res<Report> {
    if opts.precision.is_some() && host.is_none() {
        return Err("--precision needs --host: it is the precision of the contamination percent".into());
    }
    if opts.depleted_out.is_some()
        && (host.is_none() || opts.every != 1 || opts.max_reads.is_some() || opts.precision.is_some())
    {
        return Err("--depleted-out needs --host and --every 1, without --max-reads or --precision: it writes every non-host read".into());
    }
    let mut depleted = match &opts.depleted_out {
        Some(prefix) => Some(Depleted::start(prefix, r2.is_some())?),
        None => None,
    };
    let host_group = host.and_then(|h| idx.group_of[h]);
    let mut dump = match &opts.reads_out {
        Some(p) => {
            let f = std::fs::File::create(p).map_err(|e| format!("{p}: {e}"))?;
            let mut w = std::io::BufWriter::new(f);
            writeln!(w, "read_id\tcategory\tunique_hits\tshared_hits\tminimizers\tsequence")?;
            Some(w)
        }
        None => None,
    };
    let n_org = idx.names.len();
    let mut reads = vec![0u64; n_org];
    let mut hits = vec![0u64; n_org];
    let mut distinct = vec![0u64; n_org];
    let (mut own_hits, mut own_minimizers) = (vec![0u64; n_org], vec![0u64; n_org]);
    let mut special = [0u64; 6];
    let mut sampled = 0u64;
    let seen: Vec<AtomicU64> = (0..idx.len() / 64 + 1).map(|_| AtomicU64::new(0)).collect(); // one bit per index entry
    let mut ambiguous: HashMap<Vec<(usize, u32)>, u64> = HashMap::new(); // weak + conflict hit patterns, for EM
    let (mut assigned_hits, mut off_target) = (0u64, 0u64); // marker hits of assigned fragments, and those not of their organism
    let mut dups = Dups::new();
    let mut unmatched =
        Reservoir { n: opts.unmatched_n, seen: 0, rng: 0x2545F4914F6CDD1D, kept: Vec::new() };
    let (mut last_rec, mut precision_stop) = (0u64, false);

    // One thread reads R1 (a second one R2) and sends batches of sampled fragments;
    // the pool classifies a batch while the reader decompresses the next ones.
    let (tx, rx) = sync_channel::<Vec<Fragment>>(4);
    let keep_ids = dump.is_some() || opts.unmatched_out.is_some();
    let keep_raw = depleted.is_some();
    let (r1s, r2s, every, max_reads) = (r1.to_string(), r2.map(str::to_string), opts.every, opts.max_reads);
    let reader =
        std::thread::spawn(move || read_fragments(&r1s, r2s.as_deref(), every, max_reads, keep_ids, keep_raw, tx));
    let pool = rayon::ThreadPoolBuilder::new().num_threads(opts.threads).build()?;
    for batch in rx {
        let calls: Vec<Call> =
            pool.install(|| batch.par_iter().map_init(Scratch::default, |s, f| classify(idx, &seen, f, s)).collect());
        // fold in input order, so reports and --reads are the same with any number of threads
        for (f, c) in batch.iter().zip(&calls) {
            sampled += 1;
            last_rec = f.rec;
            // adapter dimers all start alike: they would form one huge "cluster"
            if let (Some(key), false) = (c.key, c.adapter) {
                dups.add(key);
            }
            let total: u32 = c.hits.iter().map(|h| h.1).sum();
            for &(g, h, d) in &c.hits {
                hits[g] += h as u64;
                distinct[g] += d as u64;
            }
            let decision = call(&c.hits, idx.genomes, &idx.group_of, c.minimizers, opts.min_frac);
            if let (Some(d), Some((raw1, raw2))) = (&mut depleted, &f.raw) {
                if !decision.is_some_and(|(g, _)| Some(g) == host || Some(g) == host_group) {
                    raw1.fastq(&f.seq, &mut d.r1);
                    if let (Some(raw2), Some(seq2)) = (raw2, &f.mate) {
                        raw2.fastq(seq2, &mut d.r2);
                    }
                }
            }
            let category = match decision {
                Some((g, support)) => {
                    reads[g] += 1;
                    own_hits[g] += c.hits.iter().find(|h| h.0 == g).map_or(0, |h| h.1) as u64;
                    own_minimizers[g] += c.minimizers as u64;
                    assigned_hits += total as u64;
                    off_target += (total - support) as u64;
                    idx.names[g].as_str()
                }
                None if !c.hits.is_empty() => {
                    // single hits are dropped: on real lanes they inflated trace organisms ~10x.
                    // EM splits among organisms only: group markers say too little.
                    let pattern: Vec<(usize, u32)> =
                        c.hits.iter().filter(|h| h.1 >= 2 && h.0 < idx.genomes).map(|h| (h.0, h.1)).collect();
                    if !pattern.is_empty() {
                        *ambiguous.entry(pattern).or_default() += 1;
                    }
                    // one organism only, but a single marker is too little evidence
                    let cat = if c.hits.len() == 1 { C_WEAK } else { C_CONFLICT };
                    special[cat] += 1;
                    SPECIAL[cat]
                }
                None => {
                    let cat = if c.shared > 0 {
                        C_SHARED
                    } else if c.adapter {
                        C_ADAPTER
                    } else if c.minimizers == 0 {
                        C_LOW // poly-G, simple repeats, or shorter than 27 bp
                    } else {
                        if opts.unmatched_out.is_some() {
                            unmatched.offer(f);
                        }
                        C_UNMATCHED
                    };
                    special[cat] += 1;
                    SPECIAL[cat]
                }
            };
            if let Some(w) = &mut dump {
                let hits: Vec<String> = c.hits.iter().map(|h| format!("{}:{}", idx.names[h.0], h.1)).collect();
                writeln!(
                    w,
                    "{}\t{category}\t{}\t{}\t{}\t{}",
                    String::from_utf8_lossy(&f.id),
                    hits.join(","),
                    c.shared,
                    c.minimizers,
                    String::from_utf8_lossy(&f.seq)
                )?;
            }
        }
        if let Some(d) = &mut depleted {
            d.flush_batch()?;
        }
        if let (Some(p), Some(h)) = (opts.precision, host) {
            let contam: u64 = (0..n_org).filter(|&g| g != h && Some(g) != host_group).map(|g| reads[g]).sum();
            let (lo, hi) = wilson(contam as f64 / sampled as f64, sampled as f64 / dups.design_effect());
            if sampled >= MIN_FOR_PRECISION && (hi - lo) / 2.0 <= p {
                precision_stop = true;
                break; // drops the receiver: the reader stops at its next batch
            }
        }
    }
    let (n, reader_stopped) = reader.join().map_err(|_| "R1 reader thread panicked")??;
    if let Some(d) = depleted {
        d.finish()?;
    }
    if let Some(w) = &mut dump {
        w.flush()?;
    }
    if let Some(p) = &opts.unmatched_out {
        unmatched.write(p)?;
    }

    // per-hit error rate, from the fragments the strict rule assigned
    let eps = (off_target as f64 / assigned_hits.max(1) as f64).clamp(1e-4, 0.5);
    let est = em(&reads, &ambiguous, eps);
    let stop = if precision_stop {
        Some("precision")
    } else if reader_stopped {
        Some("max_reads")
    } else {
        None
    };
    Ok(Report {
        host,
        host_group,
        reads,
        hits,
        distinct,
        own_hits,
        own_minimizers,
        special,
        sampled,
        total: if stop.is_some() { last_rec } else { n },
        stop,
        design: dups.design_effect(),
        eps,
        est,
    })
}

/// `bycatch scan`: one sample, report on stdout.
pub fn scan(idx: &Index, host: Option<&str>, opts: &Opts, r1: &str, r2: Option<&str>, mqc: Option<(&str, &str)>) -> Res<()> {
    let report = scan_sample(idx, opts, host_id(idx, host)?, r1, r2)?;
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    report.write_tsv(idx, &mut out)?;
    out.flush()?;
    if let Some((prefix, name)) = mqc {
        crate::run::write_mqc(prefix, idx, &[(name.to_string(), &report)])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{diversity, em, wilson, Dups};
    use std::collections::HashMap;

    #[test]
    fn wilson_values() {
        // 0 of 100: upper bound z²/(n + z²) = 3.6995 %
        let (lo, hi) = wilson(0.0, 100.0);
        assert_eq!(lo, 0.0);
        assert!((hi - 3.6995).abs() < 1e-3);
        // 50 of 100: 40.4 to 59.6 %
        let (lo, hi) = wilson(0.5, 100.0);
        assert!((lo - 40.38).abs() < 0.01 && (hi - 59.62).abs() < 0.01);
        assert!((wilson(1.0, 100.0).1 - 100.0).abs() < 1e-9);
    }

    #[test]
    fn microsatellite_starts_are_not_keys() {
        use super::informative_start;
        assert!(!informative_start(&b"GT".repeat(75)));
        assert!(!informative_start(&b"CA".repeat(75)));
        assert!(!informative_start(&[b'G'; 150]));
        assert!(!informative_start(b"ACGT"));
        assert!(informative_start(b"GACTATCCCCTGTCTCTTATACACATCTCCGA"));
    }

    #[test]
    fn duplicate_design_effect() {
        let mut d = Dups::new();
        (0..1000u64).for_each(|k| d.add(crate::kmer::mix(k)));
        assert_eq!(d.design_effect(), 1.0);
        // every fragment twice: D = 2, so the CI uses half the fragments
        (0..1000u64).for_each(|k| d.add(crate::kmer::mix(k)));
        assert!((d.design_effect() - 2.0).abs() < 1e-9);
        // past MAX keys a hash fraction is kept, and D stays about the same
        let mut d = Dups::new();
        for k in 0..6 * Dups::MAX as u64 {
            d.add(crate::kmer::mix(k / 3)); // clusters of 3
        }
        assert!(d.shift > 0 && (d.design_effect() - 3.0).abs() < 0.01, "{}", d.design_effect());
    }

    #[test]
    fn min_frac_scales_with_minimizers() {
        use super::call;
        let one = |best: u32, total: u32, minz: usize, f: f64| {
            let hits = if best == total { vec![(0, best, 0)] } else { vec![(0, best, 0), (1, total - best, 0)] };
            call(&hits, 2, &[None, None], minz, f).is_some()
        };
        assert!(one(2, 2, 40, 0.0));
        assert!(!one(2, 2, 40, 0.1), "2 of 40 minimizers is below 10 %");
        assert!(one(4, 4, 40, 0.1));
        assert!(!one(1, 1, 5, 0.0) && !one(8, 10, 40, 0.0), "2-hit and 90 % rules still apply");
    }

    #[test]
    fn groups_take_reads_their_members_share() {
        use super::call;
        // organisms 0 and 1 in group 3, organism 2 alone
        let group_of = [Some(3), Some(3), None];
        let c = |hits: &[(usize, u32, u32)]| call(hits, 3, &group_of, 20, 0.0);
        assert_eq!(c(&[(0, 5, 0), (3, 10, 0)]), Some((0, 15)), "group markers back their member");
        assert_eq!(c(&[(0, 3, 0), (1, 3, 0), (3, 10, 0)]), Some((3, 16)), "both members: the group");
        assert_eq!(c(&[(3, 6, 0)]), Some((3, 6)), "group markers only");
        assert_eq!(c(&[(2, 6, 0), (3, 6, 0)]), None, "group against an outsider: conflict");
        assert_eq!(c(&[(3, 1, 0)]), None, "one hit is too little");
    }

    #[test]
    fn em_splits_ambiguous_by_abundance() {
        // 900 reads of 0, 100 of 1; 100 conflicts hit both equally, 10 reads only 2 can explain
        let amb = HashMap::from([(vec![(0, 3), (1, 3)], 100), (vec![(2, 2)], 10)]);
        let est = em(&[900, 100, 0], &amb, 0.01);
        assert!((est.iter().sum::<f64>() - 1110.0).abs() < 1e-6);
        assert!((est[0] - 990.0).abs() < 0.5 && (est[1] - 110.0).abs() < 0.5);
        assert!((est[2] - 10.0).abs() < 1e-9, "only one candidate: all go to it");
        // more hits of 1 than of 0 outweigh 0's abundance
        let amb = HashMap::from([(vec![(0, 2), (1, 4)], 10)]);
        assert!(em(&[900, 100, 0], &amb, 0.01)[1] > 109.9);
        assert_eq!(em(&[5, 6], &HashMap::new(), 0.01), vec![5.0, 6.0]);
    }

    #[test]
    fn diversity_values() {
        assert_eq!(diversity(0, 0, 100), None);
        // 1000 hits all on different markers of a large genome: as expected
        assert!((diversity(1000, 1000, 1_000_000_000).unwrap() - 1.0).abs() < 1e-3);
        // a tiny genome saturates: 661 markers, 1.4 M hits, all 661 seen
        assert!((diversity(1_400_000, 661, 661).unwrap() - 1.0).abs() < 1e-6);
        // 16S-like pile-up: 2,500 hits on 30 markers
        assert!(diversity(2_500, 30, 569_000).unwrap() < 0.02);
    }
}
