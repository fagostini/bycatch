# bycatch compared with published contamination screening tools

As of 2026-11-18.

## Summary

bycatch is closest to Kraken2 and KrakenUniq. The other main approaches are Bracken
(abundance re-estimation), sylph (sketch-based profiling) and FCS-GX (assembly screening).

bycatch keeps canonical 27-mer window minimizers (w = 15) in a sorted, memory-mapped index.
A k-mer found in two or more genomes is tagged `shared` and never votes. A read (or pair)
is assigned to an organism when it has at least 2 marker hits and that organism holds at
least 90 % of them. A `diversity` score flags organisms whose hits pile onto a few markers.

Kraken2 and sourmash were run on bycatch's data on 2026-11-20/21 (see "Measured head-to-head"
below). Three tools written in Rust, deacon, Kun-peng and sylph, were run on the same data on
2026-11-23 (see "Measured head-to-head: Rust tools"). Bracken, centrifuger and mash are
installed but not yet run.

## Measured head-to-head: Kraken2 and sourmash (2026-11-21)

Same machine as PERFORMANCE.md (16 cores, 61 GB RAM). All three tools used the same genomes:
the 10-genome set (`genomes10.tsv`, ~9 Gb) for the simulated reads, the 23-genome set
(`genomes_all.tsv`, 31.4 Gb) for the real lane. Work files are in `/var/tmp/bench_cmp`.

- **bycatch**: current `main`, `--every 1`, `-t 16`.
- **Kraken2 2.17.2** (bioconda): default k = 35, ℓ = 31 and masking, `--threads 16`.
  - Taxonomy is flat: each genome is one species directly under the root, so a shared
    minimizer goes to the root, like bycatch's `shared`.
  - Run at `--confidence 0` (default) and `0.1`.
  - **The 2.17.1 install in `~/kraken2` is broken.** A database it builds fails to find most
    of its own minimizers: phiX against a phiX-only database matched ~3 % of its k-mers, and
    exact 150 nt substrings came back unclassified. Do not use it.
- **sourmash 4.9.4**: `sketch dna -p k=31,scaled=1000` (reads with `abund`), then
  `gather --threshold-bp 0`. Its `f_unique_weighted` counts k-mers, not reads. It is shown
  as a share of the matched k-mers, because error k-mers never match.

Every database or index was read into page cache before each scan was timed.

### Build

| | bycatch | Kraken2 | sourmash |
|---|---|---|---|
| 10 genomes: time | 112 s | 20 min (`--build`) + 16 min (add + mask) | sketching not timed separately |
| 10 genomes: size, peak RAM | 3.98 GB, 7.7 GB | 12.4 GB, 12.9 GB | — |
| 23 genomes: time | 6 min 40 s | 67 min (`--build`) + 82 min (add + mask, overlapping another build) | 4 min 16 s (8 genomes at a time, overlapping a build) |
| 23 genomes: size, peak RAM | 11 GB, 13.9 GB | 37.0 GB, 37.8 GB | 159 MB, 4.8 GB |

The Kraken2 library step also writes an uncompressed copy of every genome: 30 GB for the
23 genomes.

### Simulated reads (truth known)

`short10` has 500,000 reads of 150 nt with 0.5 % errors; 10.553 % are not human. `long10` has
20,000 long reads (mean 5.0 kb); 10.600 % are not human.

| Tool | Set | Correct | Wrong organism | Not assigned | Non-human estimate | Wall time | RSS |
|---|---|---|---|---|---|---|---|
| bycatch | short10 | 99.856 % | **0** | 0.144 % | 10.541 % | 0.8 s | 4.0 GB* |
| Kraken2 c 0 | short10 | 99.295 % | 510 (0.102 %) | 0.603 % | 10.588 % | 1.9 s | 12.7 GB |
| Kraken2 c 0.1 | short10 | 98.980 % | 92 (0.018 %) | 1.002 % | 10.472 % | 1.8 s | 12.6 GB |
| sourmash | short10 | — | — | — | 10.33 % (k-mers) | 4 s sketch + 6.5 s gather | 0.3 GB |
| bycatch | long10 | 100 % | 0 | 0 | 10.600 % | 0.5 s | 4.1 GB* |
| Kraken2 c 0 | long10 | 100 % | 0 | 0 | 10.600 % | 2.1 s | 12.6 GB |
| Kraken2 c 0.1 | long10 | 99.615 % | 0 | 0.385 % | 10.565 % | 2.1 s | 12.6 GB |
| sourmash | long10 | — | — | — | 9.78 % (16 % of k-mers matched) | 3.6 s + 6.7 s | 0.3 GB |

\* bycatch RSS is mostly the memory-mapped index in shared page cache; Kraken2 loads its
table into private memory. Kraken2 with `--memory-mapping` took 1.1 s on short10.

Kraken2's wrong calls at confidence 0 are mostly human reads called dog (237) or mouse
(213). Confidence 0.1 cuts them to 92 but leaves more reads unassigned.

### Real lane (lane A of the benchmark run, cod, 2.22 M pairs, no truth)

| Tool | Cod | Non-cod in total | Largest non-cod | Time | RSS |
|---|---|---|---|---|---|
| bycatch `--every 1` | 91.32 % | **0.053 %** (CI 0.050–0.056; EM 0.087 %) | mouse 0.024 %, frog 0.009 % | 3.9 s | 11.8 GB* |
| bycatch `--every 50` | — | 0.083 % (CI 0.060–0.115) | — | 0.9 s | 11.7 GB* |
| Kraken2 c 0 | 89.91 % | **1.45 %** (32,285 pairs) | mouse 0.36 %, dog 0.20 %, cow 0.16 %, human 0.16 % | 19.8 s | 36.9 GB |
| Kraken2 c 0.1 | 87.81 % | **0.020 %** (448 pairs) | dog 111, human 93, mouse 88 pairs | 22.8 s | 36.8 GB |
| sourmash | 99.70 % of matched k-mers | 0.30 % of matched k-mers | mouse 0.15 %, frog 0.04 %, human 0.03 % | 53 s sketch + 39 s gather (both overlapping a build) | 0.6 GB |

In this lane bycatch also reports 1.76 % `shared`, 0.97 % `conflict`, 1.69 % `weak`,
0.55 % `adapter`, 0.22 % `low_complexity` and 3.44 % `unmatched`. Kraken2 leaves 8.27 %
unclassified at confidence 0, 11.55 % at 0.1.

### Reading the results

- **Accuracy with the right reference**: on the simulated sets, all three tools land within
  0.25 points of the true 10.55 %. bycatch assigns no read to the wrong organism.
  Kraken2 at confidence 0 misassigns about 1 read in 1,000, almost all between the mammals.
- **Low-level contamination**: on the real lane the three tools disagree by a factor of 70:
  - Kraken2 at confidence 0 reports 1.45 % non-cod, spread over every mammal.
  - Kraken2 at confidence 0.1 reports 0.020 %.
  - bycatch reports 0.053 %.
  - sourmash reports 0.30 % of k-mers, which is not the same unit as reads.

  The simulated set suggests the Kraken2 c 0 figure is mostly cross-assignment of reads,
  but there is no truth for this lane. Checking a sample of each tool's non-cod reads with
  BLAST would settle it; that has not been done.
- **Speed**:
  - Classifying every read, bycatch is 2.4× faster than Kraken2 on short10 and 5× faster on
    the real lane. On the lane Kraken2 has to load 37 GB.
  - With `--every 50`, bycatch takes 0.9 s on the lane, 22× faster than Kraken2.
  - sourmash spends its time on hashing and gather, and is the slowest of the three on a
    real lane.
- **Footprint**:
  - The Kraken2 database is 3.1–3.4× larger than bycatch's index and must fit in private RAM.
  - Its build is about 20× slower (adding and masking the genomes included), and its library step needs a 30 GB uncompressed copy of the
    genomes.
  - sourmash is about 70× smaller than bycatch.
- **Still not measured**:
  - Divergent contaminants, i.e. a strain or species missing from the reference.
  - Bracken.
  - Kraken2 with a real NCBI taxonomy, where a shared minimizer can go to a genus or family
    instead of the root.

## Measured head-to-head: Rust tools (2026-11-23)

The Rust tools closest to bycatch, run on the same machine, genomes and reads as above.
bycatch, Kraken2 and sourmash were not run again: their rows are in the tables above, and
were timed on 2026-11-21. Work files are in `/var/tmp/bench_rust`; scripts are in
`data/bench/rust_cmp/`.

- **deacon 0.20.0** (cargo): a minimizer filter for host depletion, k = 31, w = 15. A read is
  a match with at least 2 distinct minimizer hits and at least 1 % of its minimizers
  (defaults). It only answers "matches this index or not", so it was run with two indexes:
  - the host (human for the simulated sets, cod for lane A);
  - a contaminant index: the union of all other genomes, minus every minimizer of the host
    (`deacon index union`, then `index diff`). This mirrors bycatch's `shared` for the host,
    but not between two contaminants.
  - Reads are counted as non-host when they hit the contaminant index. "Contaminant only"
    leaves out reads that also hit the host.
  - Pairs are scored as one fragment.
- **Kun-peng 0.7.5** (cargo): a Kraken2-style classifier with the database split into chunks.
  Defaults (k = 35, ℓ = 31, confidence 0, 2 hit groups), and the same flat taxonomy as for
  Kraken2. Run at `-T 0` and `0.1`.
  - **10 genomes only.** Its build writes temporary chunk files: the 10-genome build filled
    the last 66 GB of the disk on its first try. The 23-genome build would need about 200 GB,
    so Kun-peng was not run on lane A.
- **sylph 0.9.0**: defaults (c = 200, k = 31). `sylph query` gives the containment ANI and
  effective coverage of each genome. `sylph profile` adds abundances, but its k-mer
  reassignment step runs on one thread: on short10, and on lane A with only cod left to
  reassign, it ran 11 min without finishing, and was stopped.

### Build

| | deacon | Kun-peng | sylph |
|---|---|---|---|
| 10 genomes: time | 27 s (10 indexes) + 13 min 26 s union + 26 s diff | 9 min 7 s | 18 s |
| 10 genomes: size, peak RAM | 7.6 GB (host 2.5 + contaminant 5.1), 18.9 GB | 13.9 GB, 34.5 GB; ≥ 66 GB temporary disk | 294 MB, 2.3 GB |
| 23 genomes: time | 88 s + 37 min union + 4 min 24 s diff | not built | 41 s |
| 23 genomes: size, peak RAM | 23.0 GB (host 0.5 + contaminant 22.5), 54.2 GB | — | 879 MB, 7.2 GB |

deacon builds one genome's index fast (10 s for human), but `index union` runs on one thread
and takes most of the build time.

### Simulated reads

Same sets as above: short10 is 10.553 % non-human, long10 is 10.600 %.

| Tool | Set | Correct | Wrong organism | Not assigned | Non-human estimate | Wall time | RSS |
|---|---|---|---|---|---|---|---|
| deacon | short10 | — | 54 human reads hit the contaminant index (0.011 %) | 36 non-human reads hit neither or only human | 10.557 % (contaminant only: 10.418 %) | 4.0 s host + 8.6 s contaminant† | 4.7 / 9.5 GB |
| Kun-peng c 0 | short10 | 99.963 % | 11 (0.002 %) | 0.035 % | 10.550 % | 8.2 s | 8.4 GB |
| Kun-peng c 0.1 | short10 | 99.929 % | 2 | 0.071 % | 10.543 % | 8.1 s | 8.4 GB |
| sylph query | short10 | — | — | — | finds 8 of 10 genomes; no lsinapis, no phiX | 19.5 s | 0.6 GB |
| deacon | long10 | — | 5 human reads (0.025 %) | 0 | 10.625 % (contaminant only: 10.295 %) | 4.0 s + 8.3 s† | 5.0 / 9.6 GB |
| Kun-peng c 0 | long10 | 100 % | 0 | 0 | 10.600 % | 6.8 s | 8.4 GB |
| Kun-peng c 0.1 | long10 | 94.645 % | 0 | 5.355 % | 5.585 % | 6.5 s | 8.4 GB |
| sylph query | long10 | — | — | — | finds 4 of 10 genomes | 18.3 s | 0.4 GB |

† Most of deacon's time goes on loading the index. The simulated runs shared the machine with
a one-thread `index union`.

- Kun-peng makes 11 wrong calls on short10 at confidence 0, where Kraken2 made 510 with the
  same k, ℓ and taxonomy. One difference: the Kraken2 database was built with low-complexity
  masking, Kun-peng's without.
- At confidence 0.1 Kun-peng leaves 5.4 % of long10 unassigned, almost all from the small
  genomes: 2 of 89 phiX reads are assigned, 22 of 281 *Drosophila*, 15 of 99 *E. coli*.
  Kraken2 at 0.1 lost 0.385 %. This was not investigated.
- deacon's non-human estimate is close to the truth, but 639 non-human short reads also hit
  the human index. Without the contaminant index there is no way to tell them apart.
- sylph is built for samples with more coverage. Here every genome is at 0.005–0.08×
  effective coverage. On short10 mouse and dog get ANI 92 %, although the reads come from
  the reference itself.

### Real lane (lane A, cod, 2,223,665 pairs)

| Tool | Cod | Non-cod in total | Largest non-cod | Time | RSS |
|---|---|---|---|---|---|
| deacon | 96.03 % hit cod | 1.169 % hit the contaminant index; **0.057 %** contaminant only (1,277 pairs) | mouse 307, human 167, common frog 148 pairs (R1 of contaminant-only pairs) | 2.5 s host + 39 s contaminant | 1.2 / 37.8 GB |
| sylph query | ANI 99.50 %, 0.61× | no other genome detected | — | 3.8 s sketch + 5.4 s query | 1.0 GB |

For reference, bycatch reports 0.053 % non-cod (mouse 0.024 %, frog 0.009 %), and Kraken2
1.45 % at confidence 0 and 0.020 % at 0.1.

- deacon: 24,707 pairs (1.11 %) hit both cod and the contaminant index. Counted as
  contamination, they put deacon at Kraken2 c 0's level. Without them, deacon lands at
  bycatch's level, with the same top organisms. 86,985 pairs (3.91 %) hit neither index.
- sylph at default settings reports only cod. A trace of 0.05 % of a 2.7 Gb genome is far
  below the coverage sylph needs.

### Reading the Rust results

- **deacon** is the fastest to scan once its index is loaded: cod on lane A in 2.5 s. But it
  is a host filter. A per-organism report needs one index per organism, and composing
  indexes is slow (37 min, 54 GB for 23 genomes). Its 2-hit rule is like bycatch's. It has no
  90 % rule, so pairs with hits to both cod and a contaminant need a choice from the user.
- **Kun-peng** is more accurate than Kraken2 on short10 (11 wrong calls against 510), and its
  index is about the size of Kraken2's. Its build needs several times the database size in
  temporary disk. Its loss of long reads from small genomes at confidence 0.1 is unexplained.
- **sylph** is the smallest and needs no per-read work, but at the coverage of these sets it
  misses traces that bycatch, deacon and Kraken2 all see.

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

Kraken2, sourmash, deacon, Kun-peng and sylph were measured on bycatch's data (see the two
"Measured head-to-head" sections). The other tools were not.
