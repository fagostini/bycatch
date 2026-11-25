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

On 2026-11-24 all six tools were run again in a pre-registered benchmark on one shared
cod-host panel, with synthetic reads of known composition from non-indexed assemblies (see
"Fair benchmark"). Its accuracy numbers replace those of the two older sections.

## Tools at a glance

The tools compared in this document, described by design and documented features. Measured
speed, resources and accuracy are in "Fair benchmark" below. They are not used here.

### Strengths and weaknesses

| Tool | What it does | Strengths | Weaknesses |
|---|---|---|---|
| **bycatch** | Assigns each read or pair to one of the organisms in a fixed panel, without a taxonomy, and reports the share of each. | Exact canonical-minimizer matches, so no hash collisions. Explicit `shared` category instead of a tree; optional groups (e.g. a genus) collect what their members share. Host-depleted FASTQ output (`--depleted-out`). Categories for adapters, low-complexity reads and unmatched reads. `--every N` sampling and `--precision` early stop. Confidence intervals corrected for duplicates, plus an EM estimate. Run-level QC across samples (background, elevated samples, index hopping). MultiQC output. Memory-mapped index shared by concurrent scans. Genomes can be added to an index without a rebuild. | Hierarchy only through user-defined groups (one level, no taxonomy): reads shared outside a group are `shared`. Exact matching misses strains that differ from the reference. At most 1023 organisms and groups, and the index must fit in RAM. A new tool with a small user base. |
| **Kraken2** | Assigns reads to taxa by lowest common ancestor (LCA) on hashed minimizers, using the NCBI or a custom taxonomy. | The established standard, with a large user base and many prebuilt databases. LCA places shared reads at genus or family level. A `--confidence` threshold. Bracken re-estimates abundance from its output. | Needs a taxonomy to build. The compact hash table allows a small rate of collisions. The whole database is loaded into RAM, or memory-mapped. No sampling, intervals or run-level view. |
| **Kun-peng** | Kraken2-style LCA classifier, with a database split into chunks on disk. | Kraken2-compatible output and confidence threshold. Built for databases larger than RAM: chunks are loaded in turn. | Needs a taxonomy. A build writes large temporary files. Requires at least 3 threads. A young project. No sampling, intervals or run-level view. |
| **deacon** | Minimizer filter: a read either matches the index or not, and is kept or removed. | Writes host-depleted FASTQ/FASTA, compressed, and streams from stdin and to stdout. Absolute (`-a`) and relative (`-r`) match thresholds. Index set operations (`union`, `diff`). Designed for host depletion, with a prebuilt human pangenome index. | Does not name organisms. Separating host, contaminant and unknown takes two indexes and two passes. No abundance report beyond a JSON summary. |
| **sylph** | Sample-level profiler: containment ANI and abundance from FracMinHash sketches. | Very small databases. Scales to many thousands of genomes. Reports ANI and effective coverage. Statistical correction for low coverage. | No read-level output. Sketching samples only part of the k-mers, so it has a detection floor at low coverage. Abundance is relative to the genomes detected. |
| **sourmash** | Sample-level containment of FracMinHash sketches (`gather`). | Very small sketches. Mature, well documented, with a broad set of commands. Works with large public sketch collections. | No read-level output. Results are shares of k-mers, not reads. A minimum-overlap threshold (`--threshold-bp`) sets the detection floor. |

Hyper-Gen was considered and left out. It sketches whole genomes to estimate
genome-to-genome ANI, like Mash or skani, and has no read classification or containment
mode.

### Feature grid

Features some of the tools have and others lack. ✓ present, ✗ absent, NA not applicable
(sylph and sourmash work at sample level, not per read).

| Feature | bycatch | Kraken2 | Kun-peng | deacon | sylph | sourmash |
|---|---|---|---|---|---|---|
| Names an organism for each read | ✓ | ✓ | ✓ | ✗ | NA | NA |
| Per-read output | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ |
| Per-read threshold relative to read length | ✓ | ✓ | ✓ | ✓ | NA | NA |
| Places reads shared between genomes on a taxonomy (LCA) | ✗ ¹ | ✓ | ✓ | ✗ | NA | NA |
| Builds without a taxonomy | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ |
| Writes host-depleted reads | ✓ | ✗ ² | ✗ ² | ✓ | ✗ | ✗ |
| Adds genomes without a full rebuild | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ |
| Kraken-style report, for Kraken tooling | ✗ | ✓ | ✓ | ✗ | ✗ | ✗ |
| ANI estimate | ✗ | ✗ | ✗ | ✗ | ✓ | ✓ |
| Prebuilt databases | ✗ | ✓ | ✓ ³ | ✓ | ✓ | ✓ |

¹ One level only: optional user-defined groups collect k-mers their members share.
² `--classified-out` / `--unclassified-out` split reads by any hit or none, not by taxon.
KrakenTools extracts reads by taxon.
³ It converts Kraken2 databases (`kun_peng hashshard`).

### Gaps in existing tools, and what bycatch does about them

Features useful for contamination screening that none of the five other tools provide, or
that they provide only in part. Features one of them already has (a taxonomy and LCA,
host-depleted output, databases larger than RAM) are not gaps; they are in the feature grid
and the table above.

| Missing feature | Tools lacking it | bycatch | How |
|---|---|---|---|
| An answer from a sample of the reads, stopping once it is precise enough | all five: each reads or sketches every read | Addressed | `--every N` samples across the whole file. `--precision` stops once the confidence interval is narrow enough. |
| A confidence interval on the contamination fraction | all five (sylph gives intervals on ANI, not on abundance) | Addressed | 95 % Wilson intervals on every organism and on total contamination. |
| An estimate corrected for PCR and optical duplicates | all five | Addressed | Duplicates are found from read starts. The design effect widens the interval. |
| A run-level view: background contamination, outlier samples, index hopping | all five: each sample is judged alone | Addressed | `bycatch run`: one matrix per run, with `background`, `elevated` and `index_hopping` flags. |
| A contamination figure relative to a declared host | Kraken2, Kun-peng, sylph, sourmash report taxa, not host versus contaminant; deacon reports host versus not, without names | Addressed | `--host` gives the contamination percent, its interval and the organisms behind it. |
| Reasons for unassigned reads | Kraken2, Kun-peng and deacon give one "unclassified" bucket; sylph and sourmash have no read level | Addressed | Separate `shared`, `adapter`, `low_complexity`, `weak`, `conflict` and `unmatched` categories, plus `--unmatched-out` for follow-up with BLAST. |
| A warning when hits pile up on a few markers (a conserved gene such as 16S, a repeat, an amplicon) | all five (KrakenUniq, not compared here, counts distinct k-mers) | Addressed | The `diversity` score and `low_diversity` note. |
| Host depletion and contaminant identification in one pass | deacon depletes without names; Kraken2 and Kun-peng name reads but need KrakenTools to extract them | Addressed | `--depleted-out` writes the non-host reads during the same scan that reports the organisms. |
| Flagging a close relative of a panel organism that is not in the panel | all five | Not addressed | With no reference, bycatch also calls it the nearest panel organism, or nothing. |

## Fair benchmark: cod-host panel (2026-11-24)

This section replaces the two older head-to-head sections for accuracy. Those sections
simulated reads from the indexed genomes (a best case) and scored a real lane without truth.

The protocol was pre-registered in [bench/fair/PROTOCOL.md](bench/fair/PROTOCOL.md) and
committed (761fd11) before any tool ran. The commit also froze the scorer, the
output-to-category converters and the simulator. Count tables, timing files, summaries and
scripts are in `data/bench/fair/` (git-ignored). [bench/fair/aggregate.py](bench/fair/aggregate.py)
rebuilds the summaries from the count tables.

**Panel**: one 10-genome set for every tool (`data/genomes_cod10.tsv`, ~7 Gb):

- host: cod
- in panel: human, mouse, *L. sinapis*, *Drosophila*, *Arabidopsis*, yeast, *M. tuberculosis*, *E. coli*, phiX

**Reads**: ART HS25, paired, insert N(250, 80).

- Inserts shorter than the read run into the adapter (TruSeq plus poly-G).
- 2 % of reads get poly-G tails.
- Reads come from other assemblies than the indexed ones:
  - cod GCA_010882105.1
  - human CHM13
  - mouse CAST/EiJ
  - *Drosophila* A4
  - *Arabidopsis* Ler-0
  - yeast YJM789
  - Mtb CDC1551
  - *E. coli* Sakai
  - *L. sinapis* GCA_949711785.1
  - phiX: the panel copy, a best case.
- Off-panel spike-ins, in equal parts: Pacific cod, dog, common frog.

**Sets**: 2×50, 2×100 and 2×150, each in three mixtures, each with 3 seeds (27 sets).

| Mix | Pairs | Host | In panel | Off panel |
|---|---|---|---|---|
| low | 2 M | 99.9 % | 0.07 % | 0.03 % |
| mid | 0.5 M | 99 % | 0.7 % | 0.3 % |
| high | 0.5 M | 90 % | 7 % | 3 % |

Seed 1 of the 2×150 mid set is used only to pick each tool's threshold. Every number below is
the mean of seeds 2 and 3. A 10 M-pair 2×150 mid set is used only for timing.

**Tools**, each on its documented defaults plus one pre-registered threshold sweep:

| Tool | Sweep | Pick rule (seed 1) | Picked | Default |
|---|---|---|---|---|
| bycatch 4b937e8 | none | - | - | - |
| bycatch `--min-frac` (post-hoc) | 0–0.3 | macro F1 | 0.15 | 0.15 since this study |
| Kraken2 2.17.2, NCBI taxonomy | `--confidence` 0–0.5 | macro F1 | 0.2 | 0 |
| Kun-peng 0.7.5, NCBI taxonomy | `-T` 0–0.5 | macro F1 | 0.2 | 0 |
| deacon 0.20.0 (cod index + contaminant index) | `-a` 1–3 | contamination F1 | 3 | 2 |
| sourmash 4.9.4 gather | `--threshold-bp` 50000, 10000, 2000 | in-panel error | 2000 | 50000 |
| sylph 0.9.0 profile | none | - | - | - |

- Kraken2 and Kun-peng were also built on a flat taxonomy. The scores are almost the same,
  so only the NCBI rows are shown.
- Rust tools were built with `-C target-cpu=native`.

**Scoring**

- Kraken2 and Kun-peng calls inside a panel subtree count as that organism.
- A call on an ancestor counts as "higher", so genus *Gadus* is not cod.
- bycatch `shared`, `adapter` and `low_complexity` count as unassigned.
- deacon is scored only as host or contaminant, because it does not name organisms.
- sylph and sourmash are scored only at sample level.

### Host reads wrongly called contaminant (mid mix)

The share of cod reads called something else. It sets the false contamination floor.

| Tool | 2×50 | 2×100 | 2×150 |
|---|---|---|---|
| bycatch | 0.005 % | 0.022 % | 0.034 % |
| bycatch 0.15 (post-hoc) | 0.005 % | 0.013 % | 0.010 % |
| Kraken2 c0.2 | 0.034 % | 0.009 % | 0.007 % |
| Kraken2 c0 (default) | 0.13 % | 0.44 % | 0.64 % |
| Kun-peng c0.2 | 0.024 % | 0.001 % | 0.0003 % |
| Kun-peng c0 (default) | 0.43 % | 0.72 % | 0.86 % |
| deacon a3 | 0.002 % | 0.009 % | 0.011 % |
| deacon a2 (default) | 0.011 % | 0.020 % | 0.019 % |

### Host reads left unassigned (mid mix)

These reads are not called cod. This does not inflate contamination, but it lowers the
host count.

| Tool | 2×50 | 2×100 | 2×150 |
|---|---|---|---|
| bycatch | 19.5 % | 8.2 % | 7.4 % |
| bycatch 0.15 (post-hoc) | 19.5 % | 8.8 % | 8.4 % |
| Kraken2 c0.2 | 17.0 % | 14.5 % | 14.7 % |
| Kun-peng c0.2 | 13.6 % | 11.3 % | 11.9 % |
| deacon a3 | 52.8 % | 7.2 % | 5.3 % |
| deacon a2 (default) | 24.5 % | 5.7 % | 4.5 % |

### Per-organism accuracy

Macro F1 is the mean of per-read F1 over the 8 named in-panel organisms (cod and phiX
excluded). Each cell is low mix / high mix.

| Tool | 2×50 | 2×100 | 2×150 |
|---|---|---|---|
| bycatch | 0.87 / 0.91 | 0.88 / 0.97 | 0.85 / 0.97 |
| bycatch 0.15 (post-hoc) | 0.87 / 0.91 | 0.90 / 0.96 | 0.91 / 0.97 |
| Kraken2 c0.2 | 0.83 / 0.94 | 0.91 / 0.96 | 0.92 / 0.96 |
| Kraken2 c0 (default) | 0.73 / 0.94 | 0.66 / 0.94 | 0.61 / 0.93 |
| Kun-peng c0.2 | 0.85 / 0.93 | 0.94 / 0.96 | 0.94 / 0.95 |
| Kun-peng c0 (default) | 0.69 / 0.93 | 0.68 / 0.92 | 0.66 / 0.91 |

Contamination F1 scores each read as contaminant or not, so deacon can be compared too.
Each cell is low mix / mid mix.

| Tool | 2×50 | 2×100 | 2×150 |
|---|---|---|---|
| bycatch | 0.87 / 0.91 | 0.83 / 0.95 | 0.75 / 0.94 |
| bycatch 0.15 (post-hoc) | 0.87 / 0.91 | 0.87 / 0.95 | 0.89 / 0.96 |
| Kraken2 c0.2 | 0.74 / 0.92 | 0.89 / 0.95 | 0.90 / 0.96 |
| Kun-peng c0.2 | 0.78 / 0.92 | 0.94 / 0.96 | 0.95 / 0.95 |
| deacon a3 | 0.64 / 0.65 | 0.90 / 0.95 | 0.89 / 0.95 |
| deacon a2 (default) | 0.78 / 0.85 | 0.82 / 0.94 | 0.80 / 0.94 |

All four per-read tools detect all 8 organisms (≥ 10 pairs) in every mixture.

phiX recall:
- bycatch: 0.97–1.00
- Kraken2: 0.99
- Kun-peng: 0.73–0.83

### Contamination estimate

In-panel contamination estimate. Each cell is low mix (truth 0.07 %) / high mix (truth 7 %).

| Tool | 2×50 | 2×100 | 2×150 |
|---|---|---|---|
| bycatch | 0.065 / 6.0 % | 0.088 / 6.7 % | 0.103 / 6.8 % |
| bycatch 0.15 (post-hoc) | 0.065 / 6.0 % | 0.080 / 6.7 % | 0.077 / 6.7 % |
| Kraken2 c0.2 | 0.095 / 6.4 % | 0.075 / 6.6 % | 0.073 / 6.6 % |
| Kraken2 c0 (default) | 0.19 / 6.6 % | 0.50 / 7.2 % | 0.71 / 7.5 % |
| Kun-peng c0.2 | 0.084 / 6.2 % | 0.065 / 6.4 % | 0.063 / 6.3 % |
| deacon a3 | 0.036 / 3.6 % | 0.073 / 6.5 % | 0.075 / 6.6 % |
| sourmash th2000 | 0.165 / 6.1 % | 0.156 / 5.8 % | 0.169 / 5.8 % |
| sourmash th50000 (default) | 0.095 / 5.4 % | 0.114 / 5.1 % | 0.123 / 5.1 % |
| sylph profile | 0 / 2.5 % | timeout / 3.2 % | timeout / 5.6 % |

- The sourmash estimate is a share of k-mers, not reads.
- sylph's `Sequence_abundance` is normalised over the genomes it detects.
- At the sourmash default, the mid mix (truth 0.7 %) reads 0 / 0.32 / 0.50 %. The default
  finds 0–3 of 9 organisms at low and mid. th2000 finds all 9.
- sylph finds 0 organisms at mid and 2.5–5 of 9 at high. The 2×100 and 2×150 low-mix
  profiles did not finish within 30 min.

### Off-panel reads (mid mix)

| Tool | Called cod: 2×50 / 2×100 / 2×150 |
|---|---|
| bycatch | 21 / 31 / 32 % |
| bycatch 0.15 (post-hoc) | 21 / 30 / 30 % |
| Kraken2 c0.2 | 25 / 28 / 28 % |
| Kun-peng c0.2 | 25 / 28 / 27 % |
| deacon a3 | 8 / 30 / 33 % |

Pacific cod is one third of the off-panel reads. Every tool calls most Pacific cod reads
cod, and nearly all dog and frog reads unassigned. Off-panel reads hardly ever go to a
wrong panel organism, but not never:
- bycatch: 1.1–2.0 % at 2×100 and 2×150 (0.3 % at 2×50); 0.3–0.8 % at `--min-frac 0.15`.
- Kraken2 c0.2: 0.3–1.2 %.
- Kun-peng c0.2: 0.1–0.7 %. At its default it sends 3–4 % of 2×50 off-panel reads to one.

### Resources

16 cores and 61 GB RAM, the same machine as PERFORMANCE.md. Index and reads were read into
page cache before each run.

Build:

| Tool | Build wall | Peak RAM | Index on disk | Notes |
|---|---|---|---|---|
| bycatch | 1:47 | 5.9 GB | 3.2 GB | |
| Kraken2 | 16:39 + 18:15 | 10.1 GB | 9.7 GB | library + build |
| Kun-peng | 2:27 | 4.0 GB | 11.3 GB | 56.9 GB peak temp disk (flat build) |
| deacon | 31:44 | 18.9 GB | 6.0 GB | mostly `index union` |
| sylph | 0:16 | 1.7 GB | 0.23 GB | |
| sourmash | 3:17 | | 0.04 GB | |

Classify, 10 M pairs, 2×150 (median of 3 at 16 threads):

| Tool | 16 threads | 4 threads | 1 thread | Peak RSS | Pairs/s at 16 threads |
|---|---|---|---|---|---|
| bycatch | 20.0 s | 29.4 s | 78.6 s | 3.2 GB | 500 k |
| Kraken2 | 47.1 s | 49.0 s | 128.5 s | 9.7 GB | 212 k |
| Kun-peng | 79.3 s | 136.0 s | 216.2 s (3 threads) | 10.7 GB | 126 k |
| deacon, cod + contaminant passes | 8.8 + 14.8 s | 8.4 + 14.9 s | 19.5 + 27.0 s | 1.2 / 9.0 GB | 424 k |

- Kraken2 `--memory-mapping`: 47.3 s.
- deacon only splits host and contaminant. bycatch, Kraken2 and Kun-peng name organisms.
- `--min-frac` adds one comparison per read and did not change bycatch's speed (see Post-hoc).
- sylph and sourmash were timed only per set, not on the 10 M-pair set. On a 0.5 M-pair
  2×150 set:
  - sylph: sketch 0.9 s, `query` 2.9 s, `profile` 15:07. `profile` ran on one core despite
    `-t 16`, and timed out at 30 min on 2 M pairs.
  - sourmash: sketch 6.1 s + gather 13.1 s.

### Reading the fair benchmark

- **Thresholds matter more than tools.**
  - At their defaults, Kraken2 and Kun-peng call 0.1–0.9 % of cod reads contaminant. That
    is above the low mix's whole 0.07 % contamination.
  - At confidence 0.2, both are as clean as bycatch or cleaner from 2×100 up.
  - bycatch had no threshold, and its default sat near their tuned setting.
  - With `--min-frac 0.15` (post-hoc, now the default), bycatch comes close to tuned Kraken2
    and Kun-peng at 2×100 and 2×150 in the low mix: contamination F1 0.87–0.89 against
    0.89–0.95, and host false positives 0.010–0.013 % against 0.0003–0.009 %.
- **bycatch is best at 2×50.** Its host false positive rate (0.005 %) and contamination F1
  (0.87 low) are the best of all tools at 2×50.
- **Kraken2 and Kun-peng (tuned) lead at 2×100 and 2×150 in the low mix.** Their F1 is
  0.89–0.95 against bycatch's 0.75–0.83. bycatch's host false positive rate grows with read
  length (0.005 → 0.034 %): longer reads carry more chances for a stray exact 27-mer.
- **High mix: about even.** bycatch has the best macro F1 at 2×100 and 2×150 (0.97).
- **bycatch is the fastest full classifier.** It is 2.4× faster than Kraken2 and 4× faster
  than Kun-peng at 16 threads, with a third of the RAM and index size. deacon is as fast but
  only splits host and contaminant.
- **bycatch builds fastest but one.** Only sylph builds faster.
- **Close relatives.** No tool tells Pacific cod from Atlantic cod. Kraken2 and Kun-peng
  with the NCBI taxonomy did not move those reads to genus *Gadus*.
- **The sketch profilers are not a fit at these levels.** At their defaults, sylph and
  sourmash miss most contaminants below 1 %. sourmash th2000 detects them but
  overestimates low levels about 2×.

### Deviations from the protocol

- **Kraken2:** bioconda build only. The repository has no 2.17.2 tag to build from source.
- **Protocol location:** the protocol lives in `bench/fair/` (tracked), not `data/bench/fair/`.
- **Scoring:** contamination F1 was added to the scorer before the freeze commit, not after.
- **Kun-peng:** it asserts more than 2 threads, so its low-thread timing is at 3 threads, not 1.
- **Kun-peng temp disk:** the NCBI build's peak temp disk was lost to a bug in the disk
  poller. The flat build's figure (56.9 GB) is shown instead.
- **sylph timeout:** `profile` ran with a 30 min timeout that is not in the protocol. All
  2×100 and 2×150 low-mix profiles hit it and are missing.
- **Timing schedule:** all timings are from 2026-11-24, the same day, but at different
  hours:
  - bycatch at 13:00
  - Kraken2 at 14:00
  - Kun-peng at 18:30
  - deacon at 19:50
  - sylph/sourmash per-set times from 20:00 to 02:40

  The machine had no other heavy load.
- **Cold-cache runs skipped:** they need `drop_caches`, which needs sudo.
- **Output files:** bycatch reports and sylph outputs carry an appended timing line. The
  scorer and `aggregate.py` skip it.
- **deacon timing:** the timing runs write to `/dev/null`.

### Post-hoc: bycatch `--min-frac` (2026-11-25)

**Not part of the pre-registered run above.** The pre-registered numbers stay as they are.
The rows marked "bycatch 0.15 (post-hoc)" in the tables above come from this study.

bycatch's host false positive rate grew with read length. A study then added one option
to bycatch and tested it on new reads.

**Diagnosis.** Per-read output on the seed-1 sets showed what the false calls looked like.
At 2×100 and 2×150, cod reads called contaminant had:
- 2 marker hits, against ~20 minimizers;
- no cod hit;
- a median share of minimizers hit of 0.11–0.18.

True contaminant calls had a median share of 0.83–0.88. Only 1–2 % of them fell below 0.2.
At 2×50 the two groups do not separate (a share of 0.29–1.0 for both), but false calls are
already rare there.

**Change.**
- `--min-frac F`: the top organism must also hit at least F of the fragment's minimizers.
  This is the same idea as Kraken2's `--confidence`.
- Reads that fail go to `weak` or `conflict`.
- `--min-frac 0` gives byte-identical output to the frozen build. This was checked on one
  set. After this study the default was set to the pick, 0.15.

**Protocol, the same as for Kraken2.**
- Sweep 0, 0.05, 0.1, 0.15, 0.2 and 0.3 on 2×150 mid, seed 1. Pick by macro F1, ties going
  to 0. **Picked: 0.15.**
- Report on new seeds 4 and 5, simulated with the same scripts.
- Run Kraken2 2.17.2 (NCBI taxonomy, c0 and c0.2) on the same seeds 4 and 5.

Seeds 4–5 means:

| Mix, reads | Metric | bycatch f0 | bycatch f0.15 | Kraken2 c0.2 |
|---|---|---|---|---|
| low, 2×150 | host false positives | 0.035 % | 0.011 % | 0.008 % |
| | contamination F1 | 0.75 | 0.88 | 0.90 |
| | macro F1 | 0.86 | 0.91 | 0.91 |
| | in-panel estimate (truth 0.07 %) | 0.102 % | 0.077 % | 0.073 % |
| low, 2×100 | host false positives | 0.023 % | 0.014 % | 0.008 % |
| | contamination F1 | 0.82 | 0.87 | 0.90 |
| mid, 2×150 | contamination F1 | 0.95 | 0.96 | 0.96 |
| high, 2×150 | macro F1 | 0.967 | 0.967 | 0.958 |
| any, 2×50 | all metrics | unchanged | unchanged | |
| any | host reads left unassigned, 2×150 | 7.4 % | 8.4 % | 14.7 % |

- At 0.2, the 2×150 low-mix contamination F1 reaches 0.91 and host false positives fall to
  0.007 %. The macro F1 at the high mix drops 0.004. 0.2 was not the pick.
- Speed: 10 M pairs at 16 threads took 20.6–27.6 s per run for every setting, including the
  frozen build. The machine was noisier than on 2026-11-24, but there is no systematic cost.
  Peak memory stayed at 3.2 GB.
- With `--min-frac 0.15`, the low-mix gap to tuned Kraken2 at 2×100 and 2×150 mostly
  closes. bycatch keeps its lead in the 2×50 low mix (contamination F1 0.88 against 0.73)
  and in the 2×100 and 2×150 high mix. Kraken2 c0.2 stays ahead in the 2×50 mid and high
  mixes (macro F1 0.92 / 0.94 against 0.90 / 0.91).
- **Seeds 2–3 rows in the tables above.** The seed 2–3 sets were rebuilt from the same seeded
  scripts. bycatch at `--min-frac 0` gave count tables byte-identical to the pre-registered
  ones on all 18 sets, so the reads are the same, and so is the rebuilt seed 1 used for the
  pick. bycatch at 0.15 was then run on them.
- Unlike Kraken2's `--confidence`, this option was added after the seed 2–3 results were
  seen. That is why it is marked post-hoc and was first checked on seeds 4–5.
- Data: `data/bench/fair2/` (count tables, results, scripts; seeds 2–3 counts in
  `counts/bycatch_mf_s23/`).

## Measured head-to-head: Kraken2 and sourmash (2026-11-21)

Superseded for accuracy by "Fair benchmark" above: the reads here come from the indexed
genomes, and the real lane has no truth.

Same machine as PERFORMANCE.md (16 cores, 61 GB RAM). All three tools used the same genomes:
the 10-genome set (`genomes10.tsv`, ~9 Gb) for the simulated reads, the 23-genome set
(`genomes_all.tsv`, 31.4 Gb) for the real lane. Work files were in `/var/tmp/bench_cmp` (deleted).

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

Superseded for accuracy by "Fair benchmark" above.

The Rust tools closest to bycatch, run on the same machine, genomes and reads as above.
bycatch, Kraken2 and sourmash were not run again: their rows are in the tables above, and
were timed on 2026-11-21. Work files were in `/var/tmp/bench_rust` (deleted); scripts are in
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

For a fixed list of up to 1023 genomes and groups, bycatch is exact, lean (15.4 GB for 23 genomes,
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
- `--depleted-out` writes host-depleted FASTQ during the same scan that names the
  contaminants; deacon depletes without names, Kraken2 needs KrakenTools to extract reads.

**What published tools do better**

- Hierarchy: Kraken2/Bracken assign reads shared within a genus to the genus. bycatch does
  so only for groups named in the genome list; other shared reads stay `shared`.
- Divergence: Centrifuger, Metabuli and alignment tools still find a contaminant whose
  reference is a different strain or species. Exact 27-mers do not (caveat 6).
- Scale: sylph or ganon handle thousands of genomes. bycatch is limited to 1023 organisms
  and groups, and an index that fits in RAM.
- Statistics: sylph reports ANI and coverage. bycatch now reports 95 % Wilson intervals
  and a Bracken-like EM estimate for `weak` and `conflict` reads (see README.md), but
  `shared` reads are not redistributed.
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

Kraken2, sourmash, deacon, Kun-peng and sylph were measured on bycatch's data (see "Fair
benchmark" and the two "Measured head-to-head" sections). The other tools were not.
