# Fair benchmark protocol (pre-registered)

This file was committed before any tool ran on the benchmark sets. The scorer (`score.py`), the per-tool converters (`calls.py`) and the taxonomy mapper (`taxmap.py`) were committed with it. Anything done differently after this commit is listed under "Deviations" in COMPARISON.md.

bycatch is our own tool. This protocol exists to keep that from shaping the result.

## Frozen versions

| Tool | Version | Source |
|---|---|---|
| bycatch | commit 4b937e8 (`main`) | git worktree, `cargo build --release`, `RUSTFLAGS="-C target-cpu=native"` |
| Kraken2 | 2.17.2 | bioconda |
| deacon | 0.20.0 | `cargo install --locked`, same RUSTFLAGS |
| Kun-peng | 0.7.5 | `cargo install --locked`, same RUSTFLAGS |
| sylph | 0.9.0 | existing `/usr/local/bin/sylph` |
| sourmash | 4.9.4 | bioconda |
| ART | 2.5.8 (2016-06-06) | bioconda |

bycatch is not changed after the sets are generated.

## Reference panel

All tools index the same 10 genomes, `data/genomes_cod10.tsv`: cod (gadMor3.0), human (hg38), mouse (mm39), *L. sinapis* (ilLepSina1.1), *Drosophila* (dm6), *Arabidopsis* (TAIR10), yeast (sacCer3), *M. tuberculosis* (H37Rv), *E. coli* (K-12 MG1655), phiX. Cod is the host.

## Synthetic sets

Reads are simulated with ART (`-ss HS25`, paired, insert mean 250, sd 80) from assemblies that are **not** the indexed ones:

| Organism | Simulated from | Role |
|---|---|---|
| cod | GCA_010882105.1 (Celtic Sea individual) | host |
| human | GCF_009914755.1 (T2T-CHM13v2.0) | in panel |
| mouse | GCA_921999005.2 (CAST/EiJ v3) | in panel |
| *L. sinapis* | GCA_949711785.1 (Swedish male) | in panel |
| *Drosophila* | GCA_002300595.1 (strain A4) | in panel |
| *Arabidopsis* | GCA_965638225.1 (Ler-0) | in panel |
| yeast | GCA_000181435.1 (YJM789) | in panel |
| *M. tuberculosis* | GCF_000008585.1 (CDC1551) | in panel |
| *E. coli* | GCF_000008865.2 (O157:H7 Sakai) | in panel |
| phiX | the indexed sequence (no other exists) | in panel, best case: footnote only |
| Pacific cod | GCF_031168955.1 | off panel, close relative of the host |
| dog | `data/dog.fa.gz` | off panel |
| common frog | `data/genomes/common_frog.fa.gz` | off panel |

`mix.py` draws exact pair counts from the ART pools:

| Mix | Pairs | Host | In panel (split evenly over 9) | Off panel (split evenly over 3) |
|---|---|---|---|---|
| low | 2,000,000 | 99.9 % | 0.07 % | 0.03 % |
| mid | 500,000 | 99 % | 0.7 % | 0.3 % |
| high | 500,000 | 90 % | 7 % | 3 % |

- Read lengths: 2×50, 2×100, 2×150. Three seeds each, so 27 sets.
- Every pair gets an insert length drawn from N(250, 80). When the insert is shorter than the read, the read runs into the TruSeq adapter and then poly-G.
- 2 % of reads get a poly-G tail from a random position in their second half.
- Read names are `<organism>_<n>`, which is the truth.
- Seed 1 is used only for picking thresholds. Seeds 2 and 3 give the reported numbers.
- A separate 10 M-pair set (mid mixture, 2×150, seed 9) is used only for timing.

## Commands

All runs use 16 threads unless the row is a thread-scaling row.

- **bycatch:** `bycatch build genomes_cod10.tsv --out idx --jobs 4`, then `bycatch scan -i idx --host cod --every 1 -t 16 --reads reads.tsv R1 R2`. bycatch has no threshold option, so it runs at its fixed rule.
- **Kraken2:** library headers `kraken:taxid|<NCBI taxid>`, `kraken2-build --add-to-library` (default dust masking), `kraken2-build --build` with defaults (k 35, l 31). Classify with `kraken2 --paired --confidence C`, sweeping C over 0, 0.05, 0.1, 0.2, 0.5. The DB uses the full NCBI taxonomy. A second DB with a flat taxonomy (every genome a species under root) runs at the selected C only.
- **Kun-peng:** same library, NCBI taxonomy, `kun_peng build` defaults. Classify with `kun_peng classify -P -T C`, same C sweep. A flat-taxonomy DB runs at the selected C.
- **deacon:** defaults k 31, w 15. Host index: cod. Contaminant index: `index union` of the other 9, then `index diff` minus cod. Run `deacon filter` in keep mode, paired, with `-a A -r 0.01`, sweeping A over 1, 2, 3, once against each index.
- **sylph:** `sylph sketch` on genomes (default c 200) and on reads (`-1 R1 -2 R2`). Then `sylph profile`, with a 30 min limit per set. If it runs out of time, `sylph query` is reported instead, for detection only.
- **sourmash:** `sourmash sketch dna -p k=31,scaled=1000,abund` on genomes and on both read files merged. Then `sourmash gather`, sweeping `--threshold-bp` over 50000 (default), 10000 and 2000.

## Category mapping (`calls.py`, `taxmap.py`)

- **bycatch:** the panel organism named, or `unassigned` for `shared`, `weak`, `conflict`, `adapter`, `low_complexity` and `unmatched`.
- **Kraken2 and Kun-peng:**
  - A taxid inside a panel species' subtree → that organism.
  - An ancestor of a panel species (genus *Gadus*, Eukaryota, root …) → `higher`.
  - Anything else, unclassified, or missing from the output → `unassigned`.
- **deacon:** host index only → `cod`; contaminant index only → `contaminant`; both → `both`; neither → `unassigned`.
- **Host-vs-contamination rule:**
  - `higher` and `unassigned` count as neither host nor contamination.
  - `both` counts as host.
  - For off-panel reads, any call except `cod` and `both` counts as flagged.

## Metrics (`score.py`)

**Primary**, fixed before any run:
1. **Contamination estimate error:** |estimated − true| for the in-panel contamination fraction and for the non-host fraction.
2. **Per-organism F1:** macro average over the 8 scored in-panel organisms, phiX excluded. deacon names no organism, so it gets the binary contamination F1 instead.
3. **Host false-positive rate:** host pairs called contamination, with a 95 % Wilson interval.
4. **Off-panel handling:** shares called cod, called another panel organism, or not assigned.

**Secondary:**
- host pairs lost (not called host)
- organisms detected (≥ 10 pairs called)
- phiX recall (footnote)

**Threshold selection:**
- Kraken2 and Kun-peng: the setting with the highest macro F1 on the mid 2×150 set, seed 1.
- deacon and sourmash: the highest binary contamination F1, or the lowest in-panel estimate error for sourmash, on the same set.
- Ties go to the tool's default.
- Each tool is reported at its selected setting and at its default.

**sylph and sourmash** are scored at sample level only:
- in-panel estimate = summed abundance of non-cod panel genomes
- detected organisms

## Resources

- Same machine and same day, with nothing else heavy running.
- `/usr/bin/time -v`: wall time, CPU time, peak RSS.
- Index size, plus peak temporary disk polled every 10 s.
- Index and input are warmed with `cat` before each timed run.
- Cold-cache runs only if a page-cache drop is possible without sudo; otherwise this is skipped and stated.
- Timed runs: 3 repeats on the 10 M-pair set, median reported. Build reported apart from classify.
- Thread scaling at 1, 4 and 16 threads.

## Known limits, stated in advance

- **ART's HS25 error model is HiSeq 2500, not NovaSeq.** The NovaSeq poly-G tails are added on top.
- **The adapter read-through ignores the true fragment ends.** Each read is cut at the drawn insert length on its own, so R1 and R2 are not exact reverse complements over a short insert.
- **Kraken2 is not built from source.** Upstream has no 2.17.2 tag, and the 2.17.1 source build gave broken databases on this machine before.
- **Small contigs get no reads.** ART rounds reads per contig, so contigs shorter than about 2·L/f get none. This mostly affects the scaffold-level *L. sinapis* and yeast assemblies.
- **A 10-genome panel is a screening setup.** Kraken2 and Kun-peng are built for large reference databases, and this benchmark does not test that.
