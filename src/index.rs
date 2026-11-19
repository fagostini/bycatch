//! Reference index: sorted entries `key << 8 | organism_id`, plus a table of where each
//! k-mer prefix starts.
//!
//! The table splits the k-mers by their top `prefix_bits` bits. Large indexes (2^30
//! entries or more) use 30 prefix bits: the table then takes 4 GB, but each entry only
//! needs the remaining 24 k-mer bits plus the id, so entries shrink from 8 to 4 bytes.
//! Smaller indexes keep 8-byte entries with the full k-mer and a small table that only
//! speeds up lookups.
use crate::kmer::{kmers_at, minimizers, K, W};
use crate::scan::{adapter_kmers, ADAPTER_K};
use crate::Res;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Id of a k-mer found in more than one genome. Never evidence for any organism.
pub const SHARED: u8 = 255;

const KMER_BITS: u32 = 2 * K as u32;
/// 30 prefix bits leave 24 k-mer bits + 8 id bits: a 4-byte entry.
const COMPACT_PREFIX_BITS: u32 = 30;
/// The table stores offsets as u32, plus one u64 base per block of 2^16 prefixes.
const BLOCK_BITS: u32 = 16;

/// Prefix bits and entry size for an index of about `n` entries.
fn layout(n: u64) -> (u32, usize) {
    if n >= 1 << COMPACT_PREFIX_BITS {
        // 4 GB table + 4 bytes per entry is no larger than 8 bytes per entry from here on
        (COMPACT_PREFIX_BITS, 4)
    } else {
        // about 8 entries per prefix: table is 1/16 of the entries' size
        ((63 - n.max(16).leading_zeros()) - 3, 8)
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
fn genome_run(path: &str, tmp: &str) -> Res<u64> {
    let mut reader = needletail::parse_fastx_file(path).map_err(|e| format!("{path}: {e}"))?;
    let mut v = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| format!("{path}: {e}"))?;
        let seq = rec.seq();
        match mask_adapters(&seq) {
            Some(masked) => minimizers(&masked, &mut v),
            None => minimizers(&seq, &mut v),
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

/// TSV of `name<TAB>fasta_path`; blank lines and `#` lines are skipped.
fn read_list(list: &str, text: &str) -> Res<Vec<(String, String)>> {
    let mut genomes: Vec<(String, String)> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
        let (name, path) = line
            .split_once('\t')
            .ok_or_else(|| format!("{list}: expected name<TAB>path, got '{line}'"))?;
        if genomes.iter().any(|(n, _)| n == name) {
            return Err(format!("{list}: duplicate organism name '{name}'").into());
        }
        genomes.push((name.to_string(), path.trim().to_string()));
    }
    Ok(genomes)
}

/// `list`: TSV of `name<TAB>fasta_path`. Writes `<out>.kmers` and `<out>.meta`.
/// At most `jobs` genomes are held in memory at once; the rest goes through
/// temporary files `<out>.tmpN`, so memory does not grow with the number of genomes.
pub fn build(list: &str, out: &str, jobs: usize) -> Res<()> {
    let text = fs::read_to_string(list).map_err(|e| format!("{list}: {e}"))?;
    let genomes = read_list(list, &text)?;
    let n = genomes.len();
    if n == 0 || n > SHARED as usize {
        return Err(format!("{list}: need 1 to 255 genomes, got {n}").into());
    }
    let tmp: Vec<String> = (0..n).map(|g| format!("{out}.tmp{g}")).collect();
    let result = build_runs(&genomes, &tmp, jobs).and_then(|totals| {
        let sources = run_sources(&tmp, 0)?;
        let names: Vec<&str> = genomes.iter().map(|g| g.0.as_str()).collect();
        merge(sources, &names, &totals, out)
    });
    for t in &tmp {
        let _ = fs::remove_file(t);
    }
    result
}

/// Adds genomes to an existing index without reading its genomes again. Writes a new
/// index at `out`; the old one is left as it is. The result is the same as building
/// from the old list followed by the new one: a new genome's k-mer that is a marker of
/// an old organism, or already shared, becomes shared.
pub fn add(old: &str, list: &str, out: &str, jobs: usize) -> Res<()> {
    if old == out {
        return Err("add: --out must differ from the old index, which is read while the new one is written".into());
    }
    let idx = Index::load(old)?;
    let text = fs::read_to_string(list).map_err(|e| format!("{list}: {e}"))?;
    let genomes = read_list(list, &text)?;
    if let Some((name, _)) = genomes.iter().find(|(n, _)| idx.names.contains(n)) {
        return Err(format!("{list}: '{name}' is already in {old}").into());
    }
    let n = idx.names.len() + genomes.len();
    if genomes.is_empty() || n > SHARED as usize {
        return Err(format!("{list}: need 1 or more genomes and at most 255 in all, got {n}").into());
    }
    let tmp: Vec<String> = (0..genomes.len()).map(|g| format!("{out}.tmp{g}")).collect();
    let result = build_runs(&genomes, &tmp, jobs).and_then(|new_totals| {
        let mut sources: Vec<Box<dyn Source + '_>> = vec![Box::new(OldIndex { idx: &idx, pos: 0, prefix: 0 })];
        sources.extend(run_sources(&tmp, idx.names.len())?);
        let names: Vec<&str> = idx.names.iter().chain(genomes.iter().map(|g| &g.0)).map(String::as_str).collect();
        let totals: Vec<u64> = idx.totals.iter().chain(&new_totals).copied().collect();
        merge(sources, &names, &totals, out)
    });
    for t in &tmp {
        let _ = fs::remove_file(t);
    }
    result
}

/// Sorted (k-mer, organism id) pairs, each k-mer once.
trait Source {
    fn next(&mut self) -> Res<Option<(u64, u8)>>;
}

/// A per-genome run file: every k-mer belongs to one organism.
struct Run {
    r: BufReader<fs::File>,
    id: u8,
}

impl Source for Run {
    fn next(&mut self) -> Res<Option<(u64, u8)>> {
        Ok(next_u64(&mut self.r)?.map(|k| (k, self.id)))
    }
}

fn run_sources(tmp: &[String], first_id: usize) -> Res<Vec<Box<dyn Source + 'static>>> {
    tmp.iter()
        .enumerate()
        .map(|(g, t)| {
            let r = BufReader::with_capacity(1 << 20, fs::File::open(t).map_err(|e| format!("{t}: {e}"))?);
            Ok(Box::new(Run { r, id: (first_id + g) as u8 }) as Box<dyn Source>)
        })
        .collect()
}

/// The entries of an existing index, in order, with their full k-mers.
struct OldIndex<'a> {
    idx: &'a Index,
    pos: usize,
    prefix: usize, // prefix of the entry at `pos`, for 4-byte entries
}

impl Source for OldIndex<'_> {
    fn next(&mut self) -> Res<Option<(u64, u8)>> {
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
            (self.prefix as u64) << (KMER_BITS - idx.prefix_bits) | e >> 8
        } else {
            e >> 8
        };
        self.pos += 1;
        Ok(Some((kmer, e as u8)))
    }
}

fn build_runs(genomes: &[(String, String)], tmp: &[String], jobs: usize) -> Res<Vec<u64>> {
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
                        done.push((g, genome_run(&genomes[g].1, &tmp[g])));
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
/// k-mer coming from two sources is in two genomes and becomes `SHARED`; one from a
/// single source keeps its id (which may already be `SHARED`).
fn merge(mut runs: Vec<Box<dyn Source + '_>>, names: &[&str], totals: &[u64], out: &str) -> Res<()> {
    let mut heap = BinaryHeap::new();
    let mut ids = vec![0u8; runs.len()];
    for (g, r) in runs.iter_mut().enumerate() {
        if let Some((k, id)) = r.next()? {
            ids[g] = id;
            heap.push(Reverse((k, g)));
        }
    }
    // shared k-mers collapse, so the final count is a little lower: close enough for the layout
    let (prefix_bits, entry_bytes) = layout(totals.iter().sum());
    let suffix_bits = KMER_BITS - prefix_bits;
    let key_mask = if entry_bytes == 4 { (1u64 << suffix_bits) - 1 } else { u64::MAX };
    let path = format!("{out}.kmers");
    let mut w = BufWriter::with_capacity(1 << 20, fs::File::create(&path).map_err(|e| format!("{path}: {e}"))?);
    let mut table = TableWriter::create(&format!("{out}.table"), prefix_bits)?;
    let mut unique = vec![0u64; names.len()];
    while let Some(Reverse((kmer, g))) = heap.pop() {
        let (mut owners, first_id) = (1, ids[g]);
        if let Some((k, id)) = runs[g].next()? {
            ids[g] = id;
            heap.push(Reverse((k, g)));
        }
        while let Some(&Reverse((k2, g2))) = heap.peek() {
            if k2 != kmer {
                break;
            }
            heap.pop();
            owners += 1;
            if let Some((k, id)) = runs[g2].next()? {
                ids[g2] = id;
                heap.push(Reverse((k, g2)));
            }
        }
        let id = if owners == 1 { first_id } else { SHARED };
        if id != SHARED {
            unique[id as usize] += 1;
        }
        let id = id as u64;
        table.advance_to(kmer >> suffix_bits)?;
        let entry = (kmer & key_mask) << 8 | id;
        if entry_bytes == 4 {
            w.write_all(&(entry as u32).to_ne_bytes())?;
        } else {
            w.write_all(&entry.to_ne_bytes())?;
        }
        table.count += 1;
    }
    w.flush()?;
    table.finish()?;

    // "contam2" is the format tag from before the rename to bycatch: changing it breaks existing indexes
    let mut meta = format!("contam2\t{K}\t{W}\t{prefix_bits}\t{entry_bytes}\n");
    eprintln!("organism\tminimizers\tunique\tunique_fraction");
    for (g, name) in names.iter().enumerate() {
        meta += &format!("{name}\t{}\t{}\n", totals[g], unique[g]);
        let frac = unique[g] as f64 / totals[g].max(1) as f64;
        let warn = if frac < 0.5 { "\tWARNING: low, reduced sensitivity" } else { "" };
        eprintln!("{name}\t{}\t{}\t{frac:.3}{warn}", totals[g], unique[g]);
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

pub struct Index {
    entries: memmap2::Mmap,
    /// None for version-1 indexes: one sorted array, no prefix table
    table: Option<memmap2::Mmap>,
    prefix_bits: u32,
    entry_bytes: usize,
    /// bits of the k-mer stored in an entry
    key_mask: u64,
    len: usize,
    pub names: Vec<String>,
    /// unique marker k-mers per organism
    pub markers: Vec<u64>,
    /// minimizers per organism genome, shared ones included
    totals: Vec<u64>,
}

fn map(path: &str) -> Res<memmap2::Mmap> {
    let file = fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    // Safety: the file is opened read-only and must not be changed while we run.
    Ok(unsafe { memmap2::Mmap::map(&file)? })
}

/// First entry in `v` whose key is >= `key`, if equal: (position, id).
fn find<T: Copy + Into<u64>>(v: &[T], key: u64) -> Option<(usize, u8)> {
    let i = v.partition_point(|&e| e.into() >> 8 < key);
    (i < v.len() && v[i].into() >> 8 == key).then(|| (i, v[i].into() as u8))
}

impl Index {
    pub fn load(prefix: &str) -> Res<Index> {
        let meta_path = format!("{prefix}.meta");
        let meta = fs::read_to_string(&meta_path).map_err(|e| format!("{meta_path}: {e}"))?;
        let mut lines = meta.lines();
        let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let bad = || format!("{meta_path}: not a bycatch index for k={K}, w={W}");
        if head.len() < 3 || head[1] != K.to_string() || head[2] != W.to_string() {
            return Err(bad().into());
        }
        let (prefix_bits, entry_bytes, has_table) = match (head[0], &head[3..]) {
            ("contam1", []) => (0, 8, false),
            ("contam2", [p, b]) => (p.parse()?, b.parse()?, true),
            _ => return Err(bad().into()),
        };
        if !(entry_bytes == 8 || entry_bytes == 4 && prefix_bits == COMPACT_PREFIX_BITS) || prefix_bits > COMPACT_PREFIX_BITS {
            return Err(bad().into());
        }
        let (mut names, mut markers, mut totals) = (Vec::new(), Vec::new(), Vec::new());
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 3 {
                return Err(format!("{meta_path}: bad line '{line}'").into());
            }
            names.push(f[0].to_string());
            totals.push(f[1].parse()?);
            markers.push(f[2].parse()?);
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
        let idx = Index { entries, table, prefix_bits, entry_bytes, key_mask, len, names, markers, totals };
        if idx.table.is_some() && idx.offset(1 << prefix_bits) != len {
            return Err(format!("{prefix}.table: does not match {kmers_path}").into());
        }
        Ok(idx)
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.len
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
    pub fn lookup(&self, kmer: u64) -> Option<(usize, u8)> {
        let (lo, hi) = match self.table {
            Some(_) => {
                let p = (kmer >> (KMER_BITS - self.prefix_bits)) as usize;
                (self.offset(p), self.offset(p + 1))
            }
            None => (0, self.len),
        };
        let key = kmer & self.key_mask;
        let hit = if self.entry_bytes == 4 {
            find(&bytemuck::cast_slice::<u8, u32>(&self.entries)[lo..hi], key)
        } else {
            find(&bytemuck::cast_slice::<u8, u64>(&self.entries)[lo..hi], key)
        };
        hit.map(|(i, id)| (lo + i, id))
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

    /// `add` reads an index back through OldIndex. Checks on a real index that the decoded
    /// k-mers are sorted and look up to their own position and id:
    /// `BYCATCH_TEST_INDEX=data/all23c cargo test --release -- --ignored decode`
    #[test]
    #[ignore]
    fn decode_round_trip() {
        let prefix = std::env::var("BYCATCH_TEST_INDEX").expect("set BYCATCH_TEST_INDEX");
        let idx = Index::load(&prefix).unwrap();
        let mut src = OldIndex { idx: &idx, pos: 0, prefix: 0 };
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
