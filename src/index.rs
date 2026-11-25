//! Reference index: sorted entries `key << id_bits | organism_id`, plus a table of where
//! each k-mer prefix starts.
//!
//! The table splits the k-mers by their top `prefix_bits` bits. Large indexes use 4-byte
//! entries: the id takes as few bits as the organism count needs (7 for 95 organisms)
//! and the rest of the k-mer goes in the table, `22 + id_bits` prefix bits (2^29 slots,
//! 2 GB, for 7 id bits). Smaller indexes keep 8-byte entries with the full k-mer, 8 id
//! bits and a small table that only speeds up lookups. Over 255 ids (organisms plus
//! groups) the compact table would be 8 GB or more, so 8-byte entries are used with up to
//! 10 id bits.
//!
//! Ids: organisms first, then groups (optional third column of the genome list), then
//! the all-ones id for k-mers shared outside one group.
use crate::kmer::{kmers_at, minimizers, K, W_RANGE};
use crate::scan::{adapter_kmers, ADAPTER_K};
use crate::Res;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Id of a k-mer found in more than one genome, outside a single group. Never evidence
/// for any organism.
pub const SHARED: u16 = u16::MAX;
/// Organisms plus groups: 10 id bits next to a 54-bit k-mer fill a 64-bit entry, and
/// the all-ones id is `shared`.
pub const MAX_IDS: usize = 1023;

const KMER_BITS: u32 = 2 * K as u32;
/// k-mer bits that do not fit a 4-byte entry with no id: the least prefix bits.
const COMPACT_BASE_BITS: u32 = KMER_BITS - 32;
/// 8 id bits, the most there are, leave 24 k-mer bits in a 4-byte entry.
const COMPACT_MAX_PREFIX_BITS: u32 = COMPACT_BASE_BITS + 8;
/// The table stores offsets as u32, plus one u64 base per block of 2^16 prefixes.
const BLOCK_BITS: u32 = 16;

/// Prefix bits, entry size and id bits for an index of about `n` entries and `ids`
/// organisms plus groups.
fn layout(n: u64, ids: usize) -> (u32, usize, u32) {
    // ids 0..ids, plus the all-ones id for shared k-mers
    let id_bits = usize::BITS - ids.leading_zeros();
    let prefix_bits = COMPACT_BASE_BITS + id_bits;
    // about 8 entries per prefix: table is 1/16 of the entries' size
    let small_table = (63 - n.max(16).leading_zeros()) - 3;
    if id_bits > 8 {
        // the compact table would need 2^31 slots or more
        (small_table, 8, id_bits)
    } else if n >= 1 << prefix_bits {
        // table (4 bytes per prefix) + 4 bytes per entry is no larger than 8 bytes per entry from here on
        (prefix_bits, 4, id_bits)
    } else {
        (small_table, 8, 8)
    }
}

/// Organisms (ids 0..genomes) and groups (ids genomes..), with the group of each organism.
struct Labels {
    names: Vec<String>,
    genomes: usize,
    group_of: Vec<Option<usize>>, // per organism: id of its group
}

impl Labels {
    /// `genomes`: (name, group) in id order. Groups get ids in order of first appearance.
    fn new(genomes: &[(String, Option<String>)]) -> Res<Labels> {
        let mut names: Vec<String> = genomes.iter().map(|g| g.0.clone()).collect();
        let n = names.len();
        let mut group_of = Vec::with_capacity(n);
        for (_, group) in genomes {
            group_of.push(match group {
                None => None,
                Some(gr) => {
                    if names[..n].contains(gr) {
                        return Err(format!("group '{gr}' has the name of an organism").into());
                    }
                    Some(match names[n..].iter().position(|x| x == gr) {
                        Some(i) => n + i,
                        None => {
                            names.push(gr.clone());
                            names.len() - 1
                        }
                    })
                }
            });
        }
        if names.len() > MAX_IDS {
            return Err(format!("{n} organisms and {} groups: at most {MAX_IDS} in all", names.len() - n).into());
        }
        Ok(Labels { names, genomes: n, group_of })
    }

    /// Group id of an organism or group id, None for `shared` or an ungrouped organism.
    fn group(&self, id: u16) -> Option<usize> {
        match id as usize {
            _ if id == SHARED => None,
            g if g < self.genomes => self.group_of[g],
            g => Some(g),
        }
    }

    /// Id of a k-mer held by both `a` and `b` (organisms, groups or `shared`).
    fn combine(&self, a: u16, b: u16) -> u16 {
        match (self.group(a), self.group(b)) {
            (Some(x), Some(y)) if x == y => x as u16,
            _ => SHARED,
        }
    }
}

fn n_blocks(prefix_bits: u32) -> usize {
    ((1usize << prefix_bits) >> BLOCK_BITS) + 1
}

/// Copy of `seq` with adapter stretches set to N, or None if it has none.
/// Some assemblies (RefSeq phages, Picea abies) still end in Nextera adapter: indexed,
/// it pulls adapter read-through reads to that genome. Each matching 17-mer is widened
/// by 16 bases each side to also cover adapter remnants shorter than 17 bases.
fn mask_adapters(seq: &[u8]) -> Option<Vec<u8>> {
    let set = adapter_kmers();
    let mut out: Option<Vec<u8>> = None;
    kmers_at(seq, ADAPTER_K, |i, kmer| {
        if set.contains(&kmer) {
            let s = out.get_or_insert_with(|| seq.to_vec());
            let end = (i + 2 * ADAPTER_K - 1).min(s.len());
            s[i.saturating_sub(ADAPTER_K - 1)..end].fill(b'N');
        }
    });
    out
}

/// Sorted, deduplicated minimizers of one genome, written raw to `tmp`. Returns their number.
fn genome_run(path: &str, tmp: &str, w: usize) -> Res<u64> {
    let mut reader = needletail::parse_fastx_file(path).map_err(|e| format!("{path}: {e}"))?;
    let mut v = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| format!("{path}: {e}"))?;
        let seq = rec.seq();
        match mask_adapters(&seq) {
            Some(masked) => minimizers(&masked, w, &mut v),
            None => minimizers(&seq, w, &mut v),
        }
    }
    v.par_sort_unstable();
    v.dedup();
    fs::write(tmp, bytemuck::cast_slice::<u64, u8>(&v)).map_err(|e| format!("{tmp}: {e}"))?;
    Ok(v.len() as u64)
}

fn next_u64(r: &mut impl Read) -> Res<Option<u64>> {
    let mut b = [0u8; 8];
    match r.read_exact(&mut b) {
        Ok(()) => Ok(Some(u64::from_ne_bytes(b))),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e.into()),
    }
}

/// One line of a genome list: name, FASTA path, optional group.
type Genome = (String, String, Option<String>);

/// TSV of `name<TAB>fasta_path[<TAB>group]`; blank lines and `#` lines are skipped.
fn read_list(list: &str, text: &str) -> Res<Vec<Genome>> {
    let mut genomes: Vec<Genome> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        let (name, path, group) = match f[..] {
            [n, p] => (n, p, None),
            [n, p, g] if !g.is_empty() => (n, p, Some(g.to_string())),
            [n, p, _] => (n, p, None),
            _ => return Err(format!("{list}: expected name<TAB>path[<TAB>group], got '{line}'").into()),
        };
        if genomes.iter().any(|(n, _, _)| n == name) {
            return Err(format!("{list}: duplicate organism name '{name}'").into());
        }
        genomes.push((name.to_string(), path.to_string(), group));
    }
    Ok(genomes)
}

/// `list`: TSV of `name<TAB>fasta_path[<TAB>group]`. Writes `<out>.kmers` and `<out>.meta`.
/// At most `jobs` genomes are held in memory at once; the rest goes through
/// temporary files `<out>.tmpN`, so memory does not grow with the number of genomes.
/// `w`: minimizer window, recorded in the index; scans use the same.
pub fn build(list: &str, out: &str, jobs: usize, w: usize) -> Res<()> {
    if !W_RANGE.contains(&w) {
        return Err(format!("--window must be {} to {}, got {w}", W_RANGE.start(), W_RANGE.end()).into());
    }
    let text = fs::read_to_string(list).map_err(|e| format!("{list}: {e}"))?;
    let genomes = read_list(list, &text)?;
    let n = genomes.len();
    if n == 0 {
        return Err(format!("{list}: no genomes").into());
    }
    let labels = Labels::new(&genomes.iter().map(|g| (g.0.clone(), g.2.clone())).collect::<Vec<_>>())
        .map_err(|e| format!("{list}: {e}"))?;
    let tmp: Vec<String> = (0..n).map(|g| format!("{out}.tmp{g}")).collect();
    let result = build_runs(&genomes, &tmp, jobs, w).and_then(|totals| {
        let sources = run_sources(&tmp, 0)?;
        merge(sources, &labels, &totals, w, out)
    });
    for t in &tmp {
        let _ = fs::remove_file(t);
    }
    result
}

/// Adds genomes to an existing index without reading its genomes again. Writes a new
/// index at `out`; the old one is left as it is. The result is the same as building
/// from the old list followed by the new one: a new genome's k-mer that is a marker of
/// an old organism, or already shared, becomes shared, unless both are in one group.
pub fn add(old: &str, list: &str, out: &str, jobs: usize) -> Res<()> {
    if old == out {
        return Err("add: --out must differ from the old index, which is read while the new one is written".into());
    }
    let idx = Index::load(old)?;
    let text = fs::read_to_string(list).map_err(|e| format!("{list}: {e}"))?;
    let genomes = read_list(list, &text)?;
    if let Some((name, _, _)) = genomes.iter().find(|(n, _, _)| idx.names.contains(n)) {
        return Err(format!("{list}: '{name}' is already in {old}").into());
    }
    if genomes.is_empty() {
        return Err(format!("{list}: no genomes").into());
    }
    // old organisms keep their ids; old groups move up by the number of new organisms
    let all: Vec<(String, Option<String>)> = (0..idx.genomes)
        .map(|g| (idx.names[g].clone(), idx.group_of[g].map(|x| idx.names[x].clone())))
        .chain(genomes.iter().map(|g| (g.0.clone(), g.2.clone())))
        .collect();
    let labels = Labels::new(&all).map_err(|e| format!("{list}: {e}"))?;
    // groups of the old index that no new genome names keep their name: add them in old order
    let mut labels = labels;
    for x in idx.genomes..idx.names.len() {
        if !labels.names[labels.genomes..].contains(&idx.names[x]) {
            labels.names.push(idx.names[x].clone());
        }
    }
    let remap = old_ids(&idx, &labels);
    let tmp: Vec<String> = (0..genomes.len()).map(|g| format!("{out}.tmp{g}")).collect();
    // new genomes use the old index's window: the same k-mers must be kept
    let result = build_runs(&genomes, &tmp, jobs, idx.w).and_then(|new_totals| {
        let mut sources: Vec<Box<dyn Source + '_>> = vec![Box::new(OldIndex { idx: &idx, pos: 0, prefix: 0, remap })];
        sources.extend(run_sources(&tmp, idx.genomes)?);
        let totals: Vec<u64> = idx.totals[..idx.genomes].iter().chain(&new_totals).copied().collect();
        merge(sources, &labels, &totals, idx.w, out)
    });
    for t in &tmp {
        let _ = fs::remove_file(t);
    }
    result
}

/// Writes the index `old` again at `out` in the current layout (fewer id bits, so a
/// smaller table). Same k-mers and ids: lookups give the same results.
pub fn repack(old: &str, out: &str) -> Res<()> {
    if old == out {
        return Err("repack: --out must differ from the old index, which is read while the new one is written".into());
    }
    let idx = Index::load(old)?;
    let labels = idx.labels();
    let remap = old_ids(&idx, &labels);
    merge(vec![Box::new(OldIndex { idx: &idx, pos: 0, prefix: 0, remap })], &labels, &idx.totals[..idx.genomes], idx.w, out)
}

/// New id of each id of `idx`, by name.
fn old_ids(idx: &Index, labels: &Labels) -> Vec<u16> {
    idx.names.iter().map(|n| labels.names.iter().position(|x| x == n).unwrap() as u16).collect()
}

/// Sorted (k-mer, organism id) pairs, each k-mer once.
trait Source {
    fn next(&mut self) -> Res<Option<(u64, u16)>>;
}

/// A per-genome run file: every k-mer belongs to one organism.
struct Run {
    r: BufReader<fs::File>,
    id: u16,
}

impl Source for Run {
    fn next(&mut self) -> Res<Option<(u64, u16)>> {
        Ok(next_u64(&mut self.r)?.map(|k| (k, self.id)))
    }
}

fn run_sources(tmp: &[String], first_id: usize) -> Res<Vec<Box<dyn Source + 'static>>> {
    tmp.iter()
        .enumerate()
        .map(|(g, t)| {
            let r = BufReader::with_capacity(1 << 20, fs::File::open(t).map_err(|e| format!("{t}: {e}"))?);
            Ok(Box::new(Run { r, id: (first_id + g) as u16 }) as Box<dyn Source>)
        })
        .collect()
}

/// The entries of an existing index, in order, with their full k-mers.
struct OldIndex<'a> {
    idx: &'a Index,
    pos: usize,
    prefix: usize, // prefix of the entry at `pos`, for 4-byte entries
    remap: Vec<u16>, // new id of each old id
}

impl Source for OldIndex<'_> {
    fn next(&mut self) -> Res<Option<(u64, u16)>> {
        let idx = self.idx;
        if self.pos >= idx.len {
            return Ok(None);
        }
        let e = if idx.entry_bytes == 4 {
            bytemuck::cast_slice::<u8, u32>(&idx.entries)[self.pos] as u64
        } else {
            bytemuck::cast_slice::<u8, u64>(&idx.entries)[self.pos]
        };
        let kmer = if idx.entry_bytes == 4 {
            // the entry holds only the suffix: the prefix is the table slot it sits in
            while idx.offset(self.prefix + 1) <= self.pos {
                self.prefix += 1;
            }
            (self.prefix as u64) << (KMER_BITS - idx.prefix_bits) | e >> idx.id_bits
        } else {
            e >> idx.id_bits
        };
        self.pos += 1;
        let id = match idx.id(e) {
            SHARED => SHARED,
            id => self.remap[id as usize],
        };
        Ok(Some((kmer, id)))
    }
}

fn build_runs(genomes: &[Genome], tmp: &[String], jobs: usize, w: usize) -> Res<Vec<u64>> {
    // Plain threads, not rayon tasks: while one genome sorts, rayon work stealing
    // could start more genomes and break the memory bound.
    let next = AtomicUsize::new(0);
    let results: Vec<(usize, Res<u64>)> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..jobs.clamp(1, genomes.len()))
            .map(|_| {
                s.spawn(|| {
                    let mut done = Vec::new();
                    loop {
                        let g = next.fetch_add(1, Ordering::Relaxed);
                        if g >= genomes.len() {
                            return done;
                        }
                        eprintln!("reading {}", genomes[g].0);
                        done.push((g, genome_run(&genomes[g].1, &tmp[g], w)));
                    }
                })
            })
            .collect();
        handles.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    let mut totals = vec![0u64; genomes.len()];
    for (g, r) in results {
        totals[g] = r?;
    }
    Ok(totals)
}

/// k-way merge of sorted sources. Each source holds a k-mer at most once, so a
/// k-mer coming from two sources is in two genomes and becomes `SHARED`, or their
/// group's id if both are in one group; one from a single source keeps its id.
/// `totals`: minimizers per organism genome.
fn merge(mut runs: Vec<Box<dyn Source + '_>>, labels: &Labels, totals: &[u64], window: usize, out: &str) -> Res<()> {
    let names = &labels.names;
    let mut heap = BinaryHeap::new();
    let mut ids = vec![0u16; runs.len()];
    for (g, r) in runs.iter_mut().enumerate() {
        if let Some((k, id)) = r.next()? {
            ids[g] = id;
            heap.push(Reverse((k, g)));
        }
    }
    // shared k-mers collapse, so the final count is a little lower: close enough for the layout
    let (prefix_bits, entry_bytes, id_bits) = layout(totals.iter().sum(), names.len());
    let shared_id = (1u64 << id_bits) - 1;
    let suffix_bits = KMER_BITS - prefix_bits;
    let key_mask = if entry_bytes == 4 { (1u64 << suffix_bits) - 1 } else { u64::MAX };
    let path = format!("{out}.kmers");
    let mut w = BufWriter::with_capacity(1 << 20, fs::File::create(&path).map_err(|e| format!("{path}: {e}"))?);
    let mut table = TableWriter::create(&format!("{out}.table"), prefix_bits)?;
    let mut unique = vec![0u64; names.len()];
    while let Some(Reverse((kmer, g))) = heap.pop() {
        let mut id = ids[g];
        if let Some((k, id)) = runs[g].next()? {
            ids[g] = id;
            heap.push(Reverse((k, g)));
        }
        while let Some(&Reverse((k2, g2))) = heap.peek() {
            if k2 != kmer {
                break;
            }
            heap.pop();
            id = labels.combine(id, ids[g2]);
            if let Some((k, id)) = runs[g2].next()? {
                ids[g2] = id;
                heap.push(Reverse((k, g2)));
            }
        }
        if id != SHARED {
            unique[id as usize] += 1;
        }
        let id = if id == SHARED { shared_id } else { id as u64 };
        table.advance_to(kmer >> suffix_bits)?;
        let entry = (kmer & key_mask) << id_bits | id;
        if entry_bytes == 4 {
            w.write_all(&(entry as u32).to_ne_bytes())?;
        } else {
            w.write_all(&entry.to_ne_bytes())?;
        }
        table.count += 1;
    }
    w.flush()?;
    table.finish()?;

    // "contam2" is the format tag from before the rename to bycatch: changing it breaks
    // existing indexes. "contam3" adds the id bits and groups; it is written only when needed,
    // so indexes without groups and up to 255 organisms stay as they were.
    let v3 = id_bits > 8 || names.len() > labels.genomes;
    let mut meta = match v3 {
        true => format!("contam3\t{K}\t{window}\t{prefix_bits}\t{entry_bytes}\t{id_bits}\n"),
        false => format!("contam2\t{K}\t{window}\t{prefix_bits}\t{entry_bytes}\n"),
    };
    for x in labels.genomes..names.len() {
        eprintln!("group {}: {} unique minimizers", names[x], unique[x]);
    }
    eprintln!("organism\tminimizers\tunique\tunique_fraction");
    for (g, name) in names[..labels.genomes].iter().enumerate() {
        meta += &format!("{name}\t{}\t{}", totals[g], unique[g]);
        if v3 {
            meta += &format!("\t{}", labels.group_of[g].map_or("-", |x| names[x].as_str()));
        }
        meta += "\n";
        let frac = unique[g] as f64 / totals[g].max(1) as f64;
        let warn = if frac < 0.5 { "\tWARNING: low, reduced sensitivity" } else { "" };
        eprintln!("{name}\t{}\t{}\t{frac:.3}{warn}", totals[g], unique[g]);
    }
    for x in labels.genomes..names.len() {
        meta += &format!("{}\t0\t{}\t@group\n", names[x], unique[x]);
    }
    fs::write(format!("{out}.meta"), meta)?;
    Ok(())
}

/// Streams the prefix table. File layout: one u64 base offset per block of 2^16
/// prefixes, then a u32 offset for each prefix 0..=2^prefix_bits (the last one is the
/// number of entries). The u32 holds the offset's low 32 bits; the full offset is
/// base + (low - base as u32), wrapping, which holds while a block spans < 2^32 entries.
struct TableWriter {
    file: BufWriter<fs::File>,
    path: String,
    bases: Vec<u64>,
    /// next prefix whose start is not yet written
    next: u64,
    /// entries written so far
    count: u64,
    prefixes: u64,
}

impl TableWriter {
    fn create(path: &str, prefix_bits: u32) -> Res<TableWriter> {
        let mut file = BufWriter::with_capacity(1 << 20, fs::File::create(path).map_err(|e| format!("{path}: {e}"))?);
        // room for the block bases, filled in by finish()
        file.seek(SeekFrom::Start(8 * n_blocks(prefix_bits) as u64))?;
        Ok(TableWriter { file, path: path.into(), bases: Vec::new(), next: 0, count: 0, prefixes: 1 << prefix_bits })
    }

    /// Write the start of every prefix up to and including `prefix`.
    fn advance_to(&mut self, prefix: u64) -> Res<()> {
        while self.next <= prefix {
            if self.next.is_multiple_of(1 << BLOCK_BITS) {
                self.bases.push(self.count);
            }
            if self.count - self.bases.last().unwrap() > u32::MAX as u64 {
                return Err(format!("{}: over 2^32 entries in one table block", self.path).into());
            }
            self.file.write_all(&(self.count as u32).to_ne_bytes())?;
            self.next += 1;
        }
        Ok(())
    }

    fn finish(mut self) -> Res<()> {
        self.advance_to(self.prefixes)?;
        self.file.seek(SeekFrom::Start(0))?;
        for b in &self.bases {
            self.file.write_all(&b.to_ne_bytes())?;
        }
        self.file.flush()?;
        Ok(())
    }
}

/// `bycatch info`: layout, organisms and groups of an index, as TSV on stdout.
pub fn info(prefix: &str) -> Res<()> {
    let idx = Index::load(prefix)?;
    let meta = fs::read_to_string(format!("{prefix}.meta"))?;
    let format = meta.split('\t').next().unwrap_or("");
    let size: u64 = ["kmers", "table", "meta"]
        .iter()
        .filter_map(|e| fs::metadata(format!("{prefix}.{e}")).ok())
        .map(|m| m.len())
        .sum();
    let mut out = String::new();
    for (key, value) in [
        ("format", format.to_string()),
        ("k", K.to_string()),
        ("w", idx.w.to_string()),
        ("entry_bytes", idx.entry_bytes.to_string()),
        ("id_bits", idx.id_bits.to_string()),
        ("prefix_bits", idx.prefix_bits.to_string()),
        ("entries", idx.len.to_string()),
        ("size_bytes", size.to_string()),
        ("organisms", idx.genomes.to_string()),
        ("groups", (idx.names.len() - idx.genomes).to_string()),
    ] {
        out += &format!("# {key}\t{value}\n");
    }
    out += "organism\tminimizers\tunique\tunique_fraction\tgroup\n";
    for g in 0..idx.genomes {
        let frac = idx.markers[g] as f64 / idx.totals[g].max(1) as f64;
        let group = idx.group_of[g].map_or("-", |x| idx.names[x].as_str());
        out += &format!("{}\t{}\t{}\t{frac:.3}\t{group}\n", idx.names[g], idx.totals[g], idx.markers[g]);
    }
    if idx.names.len() > idx.genomes {
        out += "group\tunique\tmembers\n";
        for x in idx.genomes..idx.names.len() {
            let members: Vec<&str> =
                (0..idx.genomes).filter(|&g| idx.group_of[g] == Some(x)).map(|g| idx.names[g].as_str()).collect();
            out += &format!("{}\t{}\t{}\n", idx.names[x], idx.markers[x], members.join(","));
        }
    }
    print!("{out}");
    Ok(())
}

pub struct Index {
    entries: memmap2::Mmap,
    /// None for version-1 indexes: one sorted array, no prefix table
    table: Option<memmap2::Mmap>,
    prefix_bits: u32,
    entry_bytes: usize,
    /// bits of the k-mer stored in an entry
    key_mask: u64,
    /// low bits of an entry holding the organism id; all ones means shared
    id_bits: u32,
    len: usize,
    /// organisms (ids 0..genomes), then groups
    pub names: Vec<String>,
    pub genomes: usize,
    /// per organism: id of its group
    pub group_of: Vec<Option<usize>>,
    /// unique marker k-mers per organism or group
    pub markers: Vec<u64>,
    /// minimizers per organism genome, shared ones included; 0 for groups
    totals: Vec<u64>,
    /// minimizer window the index was built with
    pub w: usize,
}

fn map(path: &str) -> Res<memmap2::Mmap> {
    let file = fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    // Safety: the file is opened read-only and must not be changed while we run.
    Ok(unsafe { memmap2::Mmap::map(&file)? })
}

/// Cache hint only: never changes a result.
#[inline(always)]
fn prefetch(p: &u8) {
    #[cfg(target_arch = "x86_64")]
    // Safety: prefetch does not access memory architecturally; `p` is a valid reference anyway.
    unsafe {
        std::arch::x86_64::_mm_prefetch(p as *const u8 as *const i8, std::arch::x86_64::_MM_HINT_T0)
    };
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}

/// Position of the entry in `v` with this key, if any.
fn find<T: Copy + Into<u64>>(v: &[T], key: u64, id_bits: u32) -> Option<usize> {
    let i = v.partition_point(|&e| e.into() >> id_bits < key);
    (i < v.len() && v[i].into() >> id_bits == key).then_some(i)
}

impl Index {
    pub fn load(prefix: &str) -> Res<Index> {
        let meta_path = format!("{prefix}.meta");
        let meta = fs::read_to_string(&meta_path).map_err(|e| format!("{meta_path}: {e}"))?;
        let mut lines = meta.lines();
        let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let bad = || format!("{meta_path}: not a bycatch index for k={K}");
        if head.len() < 3 || head[1] != K.to_string() {
            return Err(bad().into());
        }
        let w: usize = head[2].parse().map_err(|_| bad())?;
        if !W_RANGE.contains(&w) {
            return Err(bad().into());
        }
        let (prefix_bits, entry_bytes, has_table, wide_bits) = match (head[0], &head[3..]) {
            ("contam1", []) => (0, 8, false, None),
            ("contam2", [p, b]) => (p.parse()?, b.parse()?, true, None),
            ("contam3", [p, b, i]) => (p.parse()?, b.parse()?, true, Some(i.parse::<u32>()?)),
            _ => return Err(bad().into()),
        };
        let compact = entry_bytes == 4 && (COMPACT_BASE_BITS + 1..=COMPACT_MAX_PREFIX_BITS).contains(&prefix_bits);
        if !(entry_bytes == 8 || compact) || prefix_bits > COMPACT_MAX_PREFIX_BITS {
            return Err(bad().into());
        }
        let id_bits = match (entry_bytes, wide_bits) {
            (4, _) => prefix_bits - COMPACT_BASE_BITS,
            (_, Some(b)) if (8..=10).contains(&b) => b,
            (_, Some(_)) => return Err(bad().into()),
            _ => 8,
        };
        let (mut names, mut markers, mut totals, mut groups) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut genomes = 0;
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            let group = match (wide_bits, &f[..]) {
                (None, [_, _, _]) => None,
                (Some(_), [_, _, _, g]) => Some(*g),
                _ => return Err(format!("{meta_path}: bad line '{line}'").into()),
            };
            if group != Some("@group") {
                if names.len() > genomes {
                    return Err(format!("{meta_path}: organism line after a group: '{line}'").into());
                }
                genomes += 1;
                groups.push(group.filter(|&g| g != "-").map(str::to_string));
            }
            names.push(f[0].to_string());
            totals.push(f[1].parse()?);
            markers.push(f[2].parse()?);
        }
        let mut group_of = Vec::with_capacity(genomes);
        for g in &groups {
            group_of.push(match g {
                None => None,
                Some(g) => Some(
                    (genomes..names.len())
                        .find(|&x| &names[x] == g)
                        .ok_or_else(|| format!("{meta_path}: unknown group '{g}'"))?,
                ),
            });
        }
        if names.len() >= 1 << id_bits {
            return Err(format!("{meta_path}: {} organisms do not fit {id_bits} id bits", names.len()).into());
        }
        let kmers_path = format!("{prefix}.kmers");
        let entries = map(&kmers_path)?;
        if entries.len() % entry_bytes != 0 {
            return Err(format!("{kmers_path}: truncated").into());
        }
        let len = entries.len() / entry_bytes;
        let table = if has_table {
            let path = format!("{prefix}.table");
            let t = map(&path)?;
            if t.len() != 8 * n_blocks(prefix_bits) + 4 * ((1 << prefix_bits) + 1) {
                return Err(format!("{path}: wrong size").into());
            }
            Some(t)
        } else {
            None
        };
        let key_mask = if entry_bytes == 4 { (1u64 << (KMER_BITS - prefix_bits)) - 1 } else { u64::MAX };
        let idx =
            Index { entries, table, prefix_bits, entry_bytes, key_mask, id_bits, len, names, genomes, group_of, markers, totals, w };
        if idx.table.is_some() && idx.offset(1 << prefix_bits) != len {
            return Err(format!("{prefix}.table: does not match {kmers_path}").into());
        }
        Ok(idx)
    }

    /// Share of an organism's minimizers that are its markers: what a read copied from
    /// its reference would hit, on average.
    pub fn unique_fraction(&self, g: usize) -> f64 {
        self.markers[g] as f64 / self.totals[g].max(1) as f64
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Organisms and groups, to write this index again.
    fn labels(&self) -> Labels {
        Labels { names: self.names.clone(), genomes: self.genomes, group_of: self.group_of.clone() }
    }

    /// Organism or group id of an entry, `SHARED` for the all-ones id.
    fn id(&self, entry: u64) -> u16 {
        let mask = (1u64 << self.id_bits) - 1;
        match entry & mask {
            m if m == mask => SHARED,
            id => id as u16,
        }
    }

    /// Position of the first entry with this k-mer prefix (needs the table).
    fn offset(&self, prefix: usize) -> usize {
        let t = self.table.as_ref().unwrap();
        let nb = 8 * n_blocks(self.prefix_bits);
        let bases: &[u64] = bytemuck::cast_slice(&t[..nb]);
        let lows: &[u32] = bytemuck::cast_slice(&t[nb..]);
        let base = bases[prefix >> BLOCK_BITS];
        (base + lows[prefix].wrapping_sub(base as u32) as u64) as usize
    }

    /// Position in the index and organism id (or `SHARED`) of a canonical k-mer.
    #[cfg(test)]
    pub fn lookup(&self, kmer: u64) -> Option<(usize, u16)> {
        let (lo, hi) = self.range(kmer);
        self.find_in(kmer, lo, hi)
    }

    /// Entries that may hold this k-mer: those with its prefix.
    pub fn range(&self, kmer: u64) -> (usize, usize) {
        match self.table {
            Some(_) => {
                let p = (kmer >> (KMER_BITS - self.prefix_bits)) as usize;
                (self.offset(p), self.offset(p + 1))
            }
            None => (0, self.len),
        }
    }

    /// Starts loading the table slot of this k-mer into cache, so that lookups of many
    /// k-mers wait for memory at the same time instead of one after another.
    pub fn prefetch_table(&self, kmer: u64) {
        if let Some(t) = &self.table {
            let p = (kmer >> (KMER_BITS - self.prefix_bits)) as usize;
            prefetch(&t[8 * n_blocks(self.prefix_bits) + 4 * p]);
        }
    }

    /// Starts loading the first entry of a `range` into cache.
    pub fn prefetch_entries(&self, lo: usize) {
        if let Some(e) = self.entries.get(lo * self.entry_bytes) {
            prefetch(e);
        }
    }

    /// `lookup` within a `range`.
    pub fn find_in(&self, kmer: u64, lo: usize, hi: usize) -> Option<(usize, u16)> {
        let key = kmer & self.key_mask;
        if self.entry_bytes == 4 {
            let v = &bytemuck::cast_slice::<u8, u32>(&self.entries)[lo..hi];
            find(v, key, self.id_bits).map(|i| (lo + i, self.id(v[i] as u64)))
        } else {
            let v = &bytemuck::cast_slice::<u8, u64>(&self.entries)[lo..hi];
            find(v, key, self.id_bits).map(|i| (lo + i, self.id(v[i])))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn masks_adapter_and_remnant() {
        let flank = b"ACGTTGCAAGGCTTACCGATGACTGGATCCATGCAATGCGTACGGATCAGTTCAGGCATAGC";
        assert!(mask_adapters(flank).is_none());
        // insert, then Nextera read-through cut 3 bases short of the 34-base adapter
        let adapter = b"CTGTCTCTTATACACATCTCCGAGCCCACGA";
        let seq = [&flank[..], adapter, flank].concat();
        let m = mask_adapters(&seq).unwrap();
        let (a, b) = (flank.len(), flank.len() + adapter.len());
        assert!(m[a..b].iter().all(|&c| c == b'N'), "adapter left: {}", String::from_utf8_lossy(&m));
        assert_eq!(&m[..a - (ADAPTER_K - 1)], &seq[..a - (ADAPTER_K - 1)]);
        assert_eq!(&m[b + ADAPTER_K - 1..], &seq[b + ADAPTER_K - 1..]);
        // reverse complement too
        let rc: Vec<u8> = seq.iter().rev().map(|&c| match c { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' }).collect();
        assert!(mask_adapters(&rc).unwrap()[a..b].iter().all(|&c| c == b'N'));
    }

    struct VecSource(std::vec::IntoIter<u64>, u16);
    impl Source for VecSource {
        fn next(&mut self) -> Res<Option<(u64, u16)>> {
            Ok(self.0.next().map(|k| (k, self.1)))
        }
    }

    /// The compact layout with few id bits, and `repack`, give the same lookups as the
    /// 8-byte layout. Fake genome totals force the layout; they only choose it.
    #[test]
    fn compact_layout_matches() {
        let dir = std::env::temp_dir().join(format!("bycatch_layout_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let p = |n: &str| dir.join(n).to_str().unwrap().to_string();
        let mut x = 7u64;
        let mut rand = || {
            x = crate::kmer::mix(x.wrapping_add(0x9e3779b97f4a7c15));
            x >> (64 - KMER_BITS)
        };
        let a: Vec<u64> = (0..20_000).map(|_| rand()).collect();
        let mut b: Vec<u64> = (0..20_000).map(|_| rand()).collect();
        b.extend(&a[..500]); // shared
        let absent: Vec<u64> = (0..1000).map(|_| rand()).collect();
        let (mut a, mut b) = (a, b);
        a.sort_unstable();
        a.dedup();
        b.sort_unstable();
        b.dedup();
        let sources = || -> Vec<Box<dyn Source>> {
            vec![Box::new(VecSource(a.clone().into_iter(), 0)), Box::new(VecSource(b.clone().into_iter(), 1))]
        };
        let names = Labels::new(&[("a".into(), None), ("b".into(), None)]).unwrap();
        merge(sources(), &names, &[10, 10], crate::kmer::W, &p("wide")).unwrap();
        merge(sources(), &names, &[1 << 24, 0], crate::kmer::W, &p("compact")).unwrap();
        let (wide, compact) = (Index::load(&p("wide")).unwrap(), Index::load(&p("compact")).unwrap());
        assert_eq!((wide.entry_bytes, compact.entry_bytes, compact.id_bits, compact.prefix_bits), (8, 4, 2, 24));
        repack(&p("compact"), &p("repacked")).unwrap();
        let repacked = Index::load(&p("repacked")).unwrap();
        let mut shared = 0;
        for &k in a.iter().chain(&b).chain(&absent) {
            let w = wide.lookup(k);
            assert_eq!(w, compact.lookup(k), "k-mer {k}");
            assert_eq!(w, repacked.lookup(k), "k-mer {k}");
            shared += (w.map(|h| h.1) == Some(SHARED)) as usize;
        }
        assert_eq!(shared, 2 * 500); // each shared k-mer is in both a and b
        assert!(absent.iter().all(|&k| wide.lookup(k).is_none()));
        fs::remove_dir_all(&dir).unwrap();
    }

    /// `add` reads an index back through OldIndex. Checks on a real index that the decoded
    /// k-mers are sorted and look up to their own position and id:
    /// `BYCATCH_TEST_INDEX=data/all23c cargo test --release -- --ignored decode`
    #[test]
    #[ignore]
    fn decode_round_trip() {
        let prefix = std::env::var("BYCATCH_TEST_INDEX").expect("set BYCATCH_TEST_INDEX");
        let idx = Index::load(&prefix).unwrap();
        let remap = (0..idx.names.len() as u16).collect();
        let mut src = OldIndex { idx: &idx, pos: 0, prefix: 0, remap };
        let (mut last, mut n) = (None, 0usize);
        while let Some((kmer, id)) = src.next().unwrap() {
            assert!(last < Some(kmer), "not sorted at {n}");
            if n % 997 == 0 {
                assert_eq!(idx.lookup(kmer), Some((n, id)), "entry {n}");
            }
            last = Some(kmer);
            n += 1;
        }
        assert_eq!(n, idx.len());
    }
}
