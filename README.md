# bycatch

Measures how much of a FASTQ file comes from each organism in a pre-selected
reference set. Built for very large files (100 GB, gzipped, short or long reads):
it reads the whole file once and classifies one read (or read pair) in N.

The reference set may mix genomes from a few Mb to 3 Gb. Sequence that two or
more genomes have in common is recognised and never counts as evidence.

Performance figures, caveats, limits and ways to scale are in
[PERFORMANCE.md](PERFORMANCE.md).

## Quick start

```bash
cargo build --release
```

Write a genome list, one organism per line: `name<TAB>fasta_path`. FASTA files may
be gzipped. Blank lines and lines starting with `#` are ignored.

```
human	data/hg38.fa.gz
fish	data/genomes/fish.fa.gz
ecoli	data/ecoli.fa
phix	data/phix.fa
```

Build the index (writes `idx.kmers`, `idx.table` and `idx.meta`):

```bash
target/release/bycatch build genomes.tsv --out idx --jobs 4
```

Scan a single-end file, a read pair, or stdin (`-`):

```bash
target/release/bycatch scan -i idx --host fish sample_R1.fastq.gz sample_R2.fastq.gz > report.tsv
```

| Option | Meaning |
|---|---|
| `-i, --index PREFIX` | index prefix given to `build --out` |
| `--host NAME` | organism the sample should be; all other assigned reads count as contamination |
| `--every N` | classify 1 read (or pair) in N, default 50; `--every 1` classifies all |
| `--max-reads M` | stop after classifying M reads (or pairs), or at the end of the file if that comes first |
| `--precision P` | stop once the contamination 95 % CI is within ± P percentage points (needs `--host`; checked every 4,096 fragments, from 10,000 on) |
| `-t, --threads T` | threads classifying reads, default 4 (plus one reading R1 and one reading R2) |
| `--reads FILE` | also write one line per classified read |
| `--unmatched-out FILE` | write a random sample of `unmatched` fragments as FASTA, for BLAST or sylph |
| `--unmatched-n N` | fragments in `--unmatched-out`, default 1000 |
| `--mqc PREFIX` | write MultiQC files `PREFIX_bycatch_mqc.tsv` and `PREFIX_bycatch_stats_mqc.tsv` |
| `--name NAME` | sample name in the MultiQC files, default from the R1 file name |
| `R1 [R2]` | FASTQ or FASTA, plain or gzipped; with R2 each pair is classified as one fragment |

`build -j/--jobs` sets how many genomes are processed at once. Each needs about
1–2 GB of RAM per Gb of genome.

`--precision` and `--max-reads` both stop early, so the sample then comes from the start
of the file only (see Reading and sampling).

### A whole run

List the run's samples, one lane per line: `sample<TAB>lane<TAB>R1[<TAB>R2]`.

```bash
target/release/bycatch run -i idx --host fish -j 8 -t 2 samples.tsv --out run1
```

`run` takes the scan options above except `--reads`, `--unmatched-out` and `--mqc`, plus
`-j/--jobs` (samples scanned at once, default 4; each holds a "seen" bitmap of one bit per
index entry, 1.4 GB for the 95-organism `count3`). It writes:

| File | Content |
|---|---|
| `run1.<sample>.<lane>.tsv` | the `scan` report of that sample and lane |
| `run1.matrix.tsv` | one row per sample and lane: sampled fragments, duplicate design effect, percent per organism and category, contamination and its CI |
| `run1.flags.tsv` | one row per sample, lane and non-host organism with ≥ 2 reads: its level, the run median, its prevalence, and flags |
| `run1_contam*_mqc.tsv` | MultiQC custom content, see below |

The flags compare each sample with the rest of the run:

| Flag | Rule | Reading |
|---|---|---|
| `background` | organism detected (≥ 2 reads) in ≥ 80 % of ≥ 5 samples, and this sample ≤ 3 × the run median | in nearly every sample at a similar level: reagent or kit contamination, or sequence missing from the host assembly, not this sample's own |
| `elevated` | ≥ 5 reads and the lower CI bound > 10 × the run median | well above the rest of the run: most likely contamination of this sample |
| `index_hopping:<sample>` | another sample on the same lane has ≥ 1 % of this organism, and this sample ≤ 2 % of that | could be reads moved from that sample by index hopping (patterned flow cells; mostly with combinatorial, not unique dual, indexes) |

A sample's own host reads hopping into other samples are not checked: there is one `--host`
for the whole run. If a sample fails (missing file, R1/R2 mismatch), the others are still
written and `run` exits with an error naming it.

MultiQC finds the `*_mqc.tsv` files when it is run on the folder holding them
(`multiqc run_folder/`; checked with MultiQC 1.35). The report gets:

| Section | From | Content |
|---|---|---|
| General Statistics | `_bycatch_stats_mqc.tsv` | Contam %, its CI upper end, Unmatched %, Adapter %, 3 decimals |
| bycatch: all fragments | `_bycatch_mqc.tsv` | stacked bars per sample, in percent: organisms (host included) and special categories |
| bycatch: contaminants | `_bycatch_contaminants_mqc.tsv` | non-host organisms only, in % of sampled fragments, so traces show next to a host at ~90 % |
| bycatch: flags | `_bycatch_flags_mqc.tsv` (`run` only) | the `elevated` and `index_hopping` rows; `background` organisms are named in the section text, as one row per sample would bury the rest |

In both bar graphs the columns are sorted by abundance (mean % over the samples). The 19
most abundant get their own colour from the Tableau 20 palette (each hue dark, then
light), and the rest are summed into one column in Tableau's dark grey, so no colour
stands for two things. In "all fragments" the host and the special categories are ranked
together with the organisms; in "contaminants" only the non-host organisms are. On the
benchmark run of [PERFORMANCE.md](PERFORMANCE.md) with the 23-genome `all23c` index this pools 9 columns in "all fragments" ("other (9)") and 2 in
"contaminants" ("other organisms (2)").

`scan --mqc` writes the first three for one sample; files from several scans with distinct
`--name`s combine into one report.

### Adding genomes

```bash
target/release/bycatch add -i idx new_genomes.tsv --out idx2
```

Reads only the new genomes and merges them into a copy of the index: their k-mers that
are markers of an old organism, or already shared, become `shared`. The result is
byte-identical to building from the old list followed by the new one (tested), and takes
the time of reading the new genomes plus one pass over the old index. The old index is
kept; `--out` must differ from it. Up to 255 organisms in all.

## How the reference index is built

The index is one sorted array of entries. Each entry holds a 27-base k-mer and
the organism it belongs to, or a tag saying it belongs to several organisms. A prefix
table says where the k-mers starting with each bit prefix begin in that array.

### 1. From sequence to k-mers

Each genome is read record by record (`needletail`; gzip detected automatically).

- **Encoding.** Each base becomes 2 bits: A = 0, C = 1, G = 2, T = 3. Lowercase
  (soft-masked) bases count as normal bases, so repeats are kept.
- **Breaks.** Any other character (`N`, IUPAC codes) ends the current stretch. No
  k-mer spans it.
- **k-mers.** Every 27-base window is a k-mer: 54 bits. k = 27 is long enough to be
  specific in genomes of several Gb, and leaves 8 bits of a `u64` for the organism id.
- **Canonical form.** The k-mer and its reverse complement are both computed while
  sliding along the sequence; the numerically smaller one is used. A read from either
  strand therefore gives the same k-mers as the genome.

### 2. Minimizers: keeping about 1 k-mer in 8

Storing every k-mer would need ~24 GB per human genome. Instead only **window
minimizers** are kept:

- Each k-mer gets a hash: the splitmix64 finalizer of its 54-bit value. The hash is
  bijective (no two k-mers share a hash) and spreads values evenly, so `AAAA…`-like
  k-mers are not favoured the way they would be by plain numeric order.
- For every run of **w = 15 consecutive k-mers** the one with the smallest hash is
  selected. Neighbouring windows usually pick the same k-mer, so on average about
  2/(w+1) ≈ 1 k-mer in 8 is kept.
- Guarantee: any stretch of k + w − 1 = **41 bases** contains at least one selected
  k-mer, and a read covering the same 41 bases selects the same one. Reads and genome
  therefore always meet on the same k-mers.
- A stretch between two `N`s that is shorter than one full window (27–40 bases) still
  contributes its single best k-mer.
- The minimum over a sliding window is kept in a monotone queue, so each base costs
  constant time.

### 3. Low-complexity filter

A selected k-mer is dropped if it contains **fewer than 8 different 3-mers** (out of
the 25 it holds). This removes homopolymers (`AAAA…`, NovaSeq poly-G) and short
tandem repeats (`ACACAC…`), which occur in every genome and in sequencing artefacts.

The same function (`kmer::minimizers`) is used for genomes and reads. If the two ever
disagreed, reads would look for k-mers the index never stored.

### 4. One sorted run per genome

For each genome, all minimizers go into one array, which is sorted (in parallel) and
deduplicated. A k-mer that occurs many times in one genome (a repeat) is now stored once.
The run is written to a temporary file `<out>.tmpN`. At most `--jobs` genomes are
in memory at once, so peak RAM depends on the largest genomes, not on how many there are.

### 5. Merge: unique markers and shared k-mers

The per-genome runs are merged in one pass (a k-way merge with a heap, reading every
run in order). For each distinct k-mer:

- found in **one** genome: entry `kmer << 8 | organism_id`. This is a **marker** of
  that organism.
- found in **two or more** genomes: one entry `kmer << 8 | 255`, the **shared** tag.

Shared entries are kept rather than deleted, so the scan can tell "sequence common to
several references" apart from "sequence in none".

This step removes the effect of conserved and contaminating sequence inside the
references themselves: a gene conserved between mouse and human, or PhiX left over in
a bacterial assembly, becomes `shared` and counts for nobody.

### 6. Entry layout and prefix table

The 54-bit k-mer is split into a **prefix** (its top P bits) and the rest. Because the
entries are sorted, all k-mers with the same prefix sit next to each other; the prefix
table stores where each group starts.

- **Large indexes** (2³⁰ entries or more, about 12 Gb of genomes and up): P = 30. The
  prefix is no longer stored in the entry, so an entry holds the remaining 24 k-mer bits
  and the 8-bit organism id: **4 bytes**. The table has 2³⁰ + 1 offsets of 4 bytes, 4 GB.
  4 GB + 4 bytes per entry beats 8 bytes per entry from 2³⁰ entries on.
- **Smaller indexes**: P is chosen so there are about 8 entries per prefix, and entries
  keep the full k-mer: `kmer << 8 | id`, **8 bytes**. The table is about 1/16 of the
  entries' size and only speeds up lookups.

Offsets are stored as 32-bit values plus one 64-bit base per block of 2¹⁶ prefixes, so
indexes with more than 2³² entries still work with a 4 GB table.

For the 23-genome index (2.77 G entries) this gives 11.1 GB of entries + 4.3 GB of table
= 15.4 GB, against 22.2 GB with 8-byte entries (−31 %). Classification is unchanged: the
same k-mers are stored with the same ids, in the same order, and reports on 16 real lane
pairs were byte-identical with both layouts.

The table also makes scans faster. Before, each lookup was a binary search over all
2.77 G entries, about 31 steps, each a likely cache miss. Now one table read narrows it to
a group of about 2–3 entries. On one paired-end lane pair (~660 Mb):

| Index | `--every 50` | `--every 1` |
|---|---|---|
| 8-byte entries, no table | 7.2 s | 196 s |
| 4-byte entries + prefix table | 3.5 s | 53 s |

Index size scales as a fixed 4.3 GB table plus 0.35 bytes per base of genome, so about
140 Gb of genomes fit in the 55 GB of free RAM on the test machine (75 Gb before).
Details in [PERFORMANCE.md](PERFORMANCE.md).

### 7. Output files

- `<out>.kmers`: the sorted entries as raw 4- or 8-byte integers, nothing else (native
  byte order, so build and scan on the same kind of machine).
- `<out>.table`: the block bases (u64), then one u32 offset per prefix.
- `<out>.meta`: text. First line `contam2	27	15	P	B` (format, k, w, prefix bits,
  bytes per entry). Then one line per organism: `name	total_minimizers	unique_markers`.
  The organism id is the line order.

Indexes from earlier versions (`contam1`: 8-byte entries, no `.table` file) still load
and give the same results, but have no prefix table, so lookups are slower.

At the end `build` prints each organism's unique fraction. Below 0.5 it prints
`WARNING: low, reduced sensitivity`: the organism shares most of its k-mers with
another genome in the list (typically two strains of one species), and few of its reads
will be assigned.

Before taking minimizers, `build` and `add` set adapter sequence in the genomes to N: every
17-mer of the scanner's adapter list (see Report below), widened by 16 bases on each side to
catch shorter remnants. Some assemblies still carry adapter at their ends (about 20 RefSeq
phage records, the *Picea abies* assembly); left in, they pulled adapter read-through reads to
that genome, and host reads with read-through became `conflict`.

## How reads are searched and classified

### 1. Reading and sampling

- The index file is memory-mapped, not loaded. Pages are read from disk on first use
  and stay in the page cache, shared by every scan running at the same time.
- Keep the index off encrypted home folders (ecryptfs): they cache every page twice,
  encrypted and decrypted, so a 48 GB index never stays in a 61 GB machine's page cache.
  A full scan of one lane took 966 s from there and 37.6 s from `/var/tmp`.
- R1 is read in full on its own thread. Only record 1, 1 + N, 1 + 2N, … is classified
  (`--every N`); the rest is decompressed and skipped. Taking every Nth read, not the
  first reads, means every part of the run (all tiles, all of the file) is represented.
- With R2, a second thread reads R2 and passes on only the sequences of the same
  record numbers. R1 and R2 must have the same number of records; otherwise the scan
  stops with an error.
- The sampled fragments go in batches of 4,096 to a pool of `--threads` classifier
  threads, while the reader goes on decompressing. Results are folded in input order,
  so the report and the `--reads` file are the same with any number of threads.
- With `--max-reads M` the scan stops right after the M-th classified fragment, so it
  reads about M × N records instead of the whole file. R2 is not read further, and the
  R1/R2 record counts are not compared. The sample then comes from the start of the file
  only (the first tiles of the flow cell), not from all of it: prefer a larger N with no
  limit when the whole run must be represented.

### 2. Looking up k-mers

For a sampled read (or pair):

1. Compute its minimizers with the same function as the build. For a pair, the
   minimizers of R1 and R2 are pooled and duplicates removed, so the pair is one fragment.
2. Look up each minimizer: its top P bits pick a slot in the prefix table, which gives
   the start and end of that prefix's group of entries (a few entries on average). A
   binary search inside the group compares the stored k-mer bits. Each lookup gives one of:
   - not found;
   - `shared`;
   - a marker of organism X, plus its position in the index.
3. Count marker hits per organism, and the number of shared hits.

### 3. Deciding the read's category

| Condition | Category |
|---|---|
| no hit, but an adapter 27-mer anywhere in the read or mate | `adapter` |
| no hit and no minimizer at all (poly-G, simple repeats, shorter than 27 bp) | `low_complexity` |
| no hit otherwise | `unmatched` |
| only `shared` hits | `shared` |
| top organism has **≥ 2** marker hits **and ≥ 90 %** of all marker hits | that organism |
| markers of one organism only, but just 1 hit | `weak` |
| anything else (markers of several organisms) | `conflict` |

Shared hits never vote. Two hits are required because one k-mer can match by chance or
through a single sequencing error. The 90 % rule tolerates an odd marker from another
organism (an error creating a k-mer that exists elsewhere) without letting a chimera be
assigned. The same rule works for long reads: they lose k-mers to errors but have
hundreds of minimizers.

The result is a fraction of **fragments**. A 3 Gb genome does not count more than a
5 Mb genome because it has more k-mers; each read is judged on its own.

### 4. Distinct markers and the `low_diversity` flag

For every marker hit, the scan also records whether that index position was hit before
(one bit per index entry). This gives, per organism, the number of **distinct** markers hit.

Reads from a real genome spread over many different markers. Reads from a conserved gene
of an organism that is not in the index (e.g. a 16S rRNA amplicon matching *E. coli*), or
from one repeated region, hit the same few markers again and again.

The `diversity` column compares the two:

    diversity = distinct / (markers × (1 − e^(−hits / markers)))

The denominator is the number of distinct markers expected if the same number of hits
fell uniformly over the organism's markers. Real genomes measured 0.59–1.0; 16S artefacts
0.005–0.16. With at least 20 hits and diversity below 0.25, the organism's row gets the
note `low_diversity`. Heavy PCR duplication can also cause it, so the note means
"check this", not "false".

### 5. Report

TSV on stdout. One row per organism, in index order, then `shared`, `weak`, `conflict`,
`adapter`, `low_complexity` and `unmatched`.

`adapter` looks for every 17-mer (not only minimizers) of reads with no other hit in a
built-in list: Illumina TruSeq read 1 and 2 and the P5 side, Nextera read 1 and 2 and the
mosaic end, the P5 and P7 flow-cell ends, MGI read 1 and 2. Exact matches only, so an insert
must leave at least 17 adapter bases in the read (insert ≤ 133 bp for 150 bp reads).

| Column | Meaning |
|---|---|
| `category` | organism name or special category |
| `reads` | sampled fragments assigned |
| `percent` | of sampled fragments |
| `marker_hits` | marker k-mer hits, counted over all sampled fragments |
| `distinct_markers` | different markers hit |
| `index_markers` | unique markers of the organism in the index |
| `diversity` | see above |
| `note` | `host` and/or `low_diversity` |
| `ci95_low`, `ci95_high` | 95 % Wilson score interval of `percent` (all rows) |
| `em_percent` | `percent` after the EM step below (organism rows only) |

The interval covers sampling noise only (`--every N`), not classification bias such as
genome reads lost to `shared`. With 0 reads its upper end, about 3/n, is the most an
organism can be present at: "phix < 0.003 %" rather than a bare 0.

**Duplicates.** PCR and optical duplicates are copies of one molecule, so they carry less
information than independent reads. Sampled fragments are grouped by a hash of the first
32 bases of the read and of its mate; with cluster sizes m, the design effect is
D = Σm² / Σm (1 without duplicates, 2 if every fragment is there twice), and the interval
uses n / D fragments. At `--every 50` two copies are rarely both sampled, so D stays near 1:
sampling already removes most of the duplication. Above 2 M distinct keys only a hash
fraction of the clusters is kept, which bounds memory and keeps D unbiased. Errors in the
first 32 bases split a cluster, so D is a slight underestimate.

**EM step.** `weak` and `conflict` fragments are split between the organisms they hit, in
proportion to abundance, as Bracken does for Kraken. Each marker hit of a fragment from
organism g is taken to be a marker of g with probability 1 − ε, else of another organism;
ε is the share of off-target hits among strictly assigned fragments (clamped to
0.0001–0.5). Only organisms with at least 2 hits in the fragment are candidates, so 1-hit
`weak` reads stay unassigned; `shared` and `unmatched` fragments are not moved. On
simulated 150 nt reads (10 genomes) it cut the largest error from −0.13 to −0.07
percentage points. On a real cod lane it raised every trace organism about 1.6×; there
is no truth to tell whether that is closer. Without the 2-hit filter it raised them up to
15×, from 1-hit reads that are most likely noise. Treat `em_percent` as a second
estimate; the strict columns are unchanged.

Trailing lines:

- `# total_fragments`: reads (or pairs) in the file, or read so far if the scan stopped early;
- `# stopped_at_max_reads`: only when `--max-reads` ended the scan before the end of the file;
- `# stopped_at_precision`: only when `--precision` ended the scan;
- `# sampled_fragments`: reads (or pairs) classified;
- `# duplicate_design_effect`: D above;
- `# em_marker_error`: the ε used by the EM step;
- `# contamination_percent`: only with `--host`; non-host assigned fragments divided by
  sampled fragments. The special categories are not counted as contamination;
- `# contamination_ci95`: only with `--host`; its 95 % Wilson interval, low and high;
- `# em_contamination_percent`: only with `--host`; the same from `em_percent`.

### 6. Per-read output (`--reads`)

One line per classified fragment:
`read_id  category  unique_hits  shared_hits  minimizers  sequence`.
`unique_hits` lists the organisms hit, such as `fish:14,mouse:1`. The sequence is the R1 sequence.
Useful to check odd calls, for example by BLAST.

`--unmatched-out FILE` writes `--unmatched-n` fragments drawn at random from all
`unmatched` ones (reservoir sampling with a fixed seed, so the same input gives the same
file). Ids are cut at the first space; pairs are written as `id/1` and `id/2`. BLAST them,
or run sylph against GTDB, to find organisms missing from the index.

## Code layout

| File | Content |
|---|---|
| `src/kmer.rs` | 2-bit encoding, canonical k-mers, hash, low-complexity filter, minimizers |
| `src/index.rs` | build and add (per-genome runs, merge, prefix table), `.meta` format, load, lookup |
| `src/scan.rs` | sampling, reader threads, parallel classification, adapters, duplicates, Wilson intervals, EM, `--precision`, report |
| `src/run.rs` | `run`: many samples, matrix, background / elevated / index-hopping flags, MultiQC files |
| `src/main.rs` | command line |
| `tests/synthetic.rs` | end-to-end test on three random genomes sharing a block: scan, threads, adapters, `--precision`, `add`, `run`, MultiQC |

```bash
cargo test --release
```
