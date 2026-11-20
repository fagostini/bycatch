# bycatch compared with published contamination screening tools

As of 2026-11-18.

## Summary

bycatch is closest to Kraken2 and KrakenUniq. The other main approaches are Bracken
(abundance re-estimation), sylph (sketch-based profiling) and FCS-GX (assembly screening).

bycatch keeps canonical 27-mer window minimizers (w = 15) in a sorted, memory-mapped index.
A k-mer found in two or more genomes is tagged `shared` and never votes. A read (or pair)
is assigned to an organism when it has at least 2 marker hits and that organism holds at
least 90 % of them. A `diversity` score flags organisms whose hits pile onto a few markers.

On this machine kraken2, bracken, centrifuger, sylph and mash are installed. None of them
has been run on bycatch's data yet.

## Per-read k-mer classifiers

These tools work like bycatch: they call each read from exact k-mer or minimizer hits. They
differ mainly in how they handle k-mers shared between genomes.

| Tool | How it differs from bycatch | Performance (approximate, published) |
|---|---|---|
| Kraken2 (Wood 2019) | Minimizers (k = 35, ℓ = 31, spaced seeds) in a compact probabilistic hash table. A k-mer shared by several genomes goes to their lowest common ancestor (LCA) in a taxonomy, not to a single `shared` tag. Reads are assigned along the root-to-leaf path, with a `--confidence` threshold. | Very fast, millions of reads per minute, multithreaded. Standard database about 50–100 GB, or capped at 8/16 GB with lower sensitivity. Hash collisions cause a small false-positive rate. Many false positives at low abundance unless `--confidence` is about 0.1 or higher. |
| KrakenUniq (Breitwieser 2018) | Kraken 1 exact k-mers, plus a HyperLogLog count of distinct k-mers per taxon. This is the same idea as bycatch's `distinct_markers` and `diversity`. | Its distinct-k-mer count is the standard published way to reject 16S and repeat artefacts. Exact database, so large (about 1 byte per base or more). Slower than Kraken2. |
| Bracken (Lu 2017) | Runs after Kraken2. A Bayesian step moves reads stuck at higher taxonomy levels down to species. | Better species-level abundance than raw Kraken output. bycatch has no such step: its `shared` and `conflict` reads are simply not counted. |
| CLARK (Ounit 2015) | Keeps only discriminative k-mers per target and deletes shared ones. | Precise, but needs a lot of RAM. bycatch keeps shared k-mers, so common sequence is not mistaken for absent sequence. |
| ganon / ganon2 (Piro 2020, 2024) | Interleaved Bloom filters plus minimizers. Shared k-mers are resolved by a re-assignment step (EM-like) instead of LCA. | Smaller indexes and faster builds than Kraken2, with similar precision. Suited to frequent RefSeq updates. |
| Taxor, KMCP, Metabuli | Taxor uses a hierarchical interleaved Bloom filter (HIBF). KMCP checks k-mer coverage across genome chunks, close to bycatch's diversity check. Metabuli also matches amino acids, for divergent organisms. | KMCP is strong on low-abundance precision. Metabuli helps when the contaminant is a different strain or species from the reference. |

## FM-index and alignment tools

These tools tolerate more sequence divergence than exact 27-mers, at a cost in speed or
index size.

| Tool | How it differs from bycatch | Performance (approximate, published) |
|---|---|---|
| Centrifuge / Centrifuger (Kim 2016, Song 2024) | BWT/FM-index with variable-length exact matches, not a fixed k. Genomes are compressed against each other. | Index about 4–10× smaller than Kraken. Somewhat more tolerant of divergence. Slower per read. Centrifuger fixes most of Centrifuge's speed and precision problems. |
| FastQ Screen, BBSplit | Align reads with Bowtie2 or BBMap against each genome. | Most tolerant of divergence. Roughly 10–100× slower than bycatch. PERFORMANCE.md already lists them as cross-checks only. |
| NCBI FCS-GX (Astashyn 2024) | Built to screen assemblies, not reads. Very large multi-kingdom database (about 470 GB). | Strongest published taxonomy-wide contamination screen, but too heavy for routine read QC. |

## Sketch and abundance tools

These tools give a presence and abundance profile per sample, with no per-read calls.

| Tool | How it differs from bycatch | Performance (approximate, published) |
|---|---|---|
| sylph (Shaw & Yu 2024) | FracMinHash (c ≈ 200). Containment ANI is corrected statistically for low coverage. | Indexes of tens to hundreds of MB. Very fast. Detects organisms at about 0.1× coverage and reports an ANI estimate. PERFORMANCE.md already names it as the best fit for hundreds of species or 50 nt reads. |
| sourmash gather, Mash screen | FracMinHash or MinHash containment against reference sketches. | Small and fast. No read-level output and no per-read evidence. |

## Where bycatch stands

For a fixed list of up to 255 genomes, bycatch is exact, lean (15.4 GB for 23 genomes,
31.4 Gb; 45.8 GB for 95 organisms, ~119 Gb) and needs no taxonomy. It falls behind on close relatives, divergent
contaminants, scale and statistics.

**What bycatch does better**

- Exact matches, so no hash false positives (unlike Kraken2).
- An explicit `shared` category, and no taxonomy needed.
- The 2-hit / 90 % rule works for both short and long reads.
- The `diversity` flag does what KrakenUniq's distinct-k-mer count does.
- A memory-mapped index shared by all scans running at once.
- `--every N` sampling over the whole file, which none of the tools above offer natively.
- Run-level QC (`bycatch run`): one matrix for all samples of a run, and flags for
  background contamination (present in most samples, as the decontam prevalence idea),
  samples well above the run median, and possible index hopping from a sample on the same
  lane. The classifiers above report each sample alone; decontam needs a separate step
  and negative controls.
- `adapter` and `low_complexity` categories, and `--unmatched-out` for the rest, instead of
  one bucket of unclassified reads.
- Confidence intervals corrected for PCR and optical duplicates (design effect), and
  `--precision` to stop reading once the contamination is known well enough.
- `bycatch add` merges new genomes into an existing index without rereading the old ones;
  Kraken2 rebuilds its whole database.
- MultiQC custom-content files from `scan --mqc` and `run`.

**What published tools do better**

- Hierarchy: Kraken2/Bracken assign reads shared within a genus to the genus. bycatch loses
  them (caveat 5 in PERFORMANCE.md).
- Divergence: Centrifuger, Metabuli and alignment tools still find a contaminant whose
  reference is a different strain or species. Exact 27-mers do not (caveat 6).
- Scale: sylph or ganon handle thousands of genomes. bycatch is limited to 255 organisms and
  an index that fits in RAM.
- Statistics: sylph reports ANI and coverage. bycatch now reports 95 % Wilson intervals
  and a Bracken-like EM estimate for `weak` and `conflict` reads (see README.md), but
  `shared` reads are still not redistributed: that needs shared-group ids in the index.
- Threads: closed. Classification now runs on `--threads` workers: one real lane pair at
  `--every 1` takes 5.1 s instead of 34.1 s with 16 threads.

**Expected speed (estimate, not measured)**

- At `--every 1`, one real lane pair takes 34.1 s on one thread and 5.1 s on 16 (6.0 s
  on 8). Above 8 threads decompression and memory bandwidth set the limit.
- At `--every 50`, bycatch takes 3.5 s per lane pair and is limited by decompression. No
  competitor is likely to be much faster there.

## Sources and caveats

The figures for the published tools are approximate. They come from memory of the papers
cited above, which were not opened or checked for this document. The figures for bycatch
come from [README.md](README.md) and [PERFORMANCE.md](PERFORMANCE.md).

No tool was benchmarked on bycatch's data. A fair head-to-head test would use the simulated
sets (`short10`, `long10`) and one real lane pair. Building a Kraken2 index of the same
31 Gb takes hours and likely more than the 61 GB of RAM on this machine, so that test may
need the 10-genome index (`real10`).
