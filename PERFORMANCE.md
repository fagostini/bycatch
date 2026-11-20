# bycatch: performance, caveats, limits and ways to scale

All numbers were measured on this workstation (16 cores, 61 GB RAM, NVMe disk)
with the current design: canonical 27-mers, window minimizers with w = 15, and a
memory-mapped, sorted index with a prefix table (4-byte entries for indexes of 2³⁰
entries or more, 8-byte entries below). Rows marked *8-byte* were measured before the
compact layout existed.

## 1. Measured performance

### Index build

| Index | Genomes | Total size | Index size | Build time | Peak RAM |
|---|---|---|---|---|---|
| `real4` (8-byte) | human, mouse, *E. coli*, PhiX | ~5.8 Gb | 4.6 GB | 55 s | 7.9 GB (old in-memory build) |
| `real10` (8-byte) | 10 genomes, 3 mammals | ~9 Gb | 7.4 GB | 90 s | 8.1 GB (11.3 GB with the old build) |
| `all23` (8-byte) | 23 genomes, 13 of them ≥ 1 Gb | 31.4 Gb | 22.2 GB | 324 s | 14.2 GB (the old build would have needed ~46 GB) |
| `all23c` (compact) | same 23 genomes | 31.4 Gb | **15.4 GB** (11.1 GB entries + 4.3 GB table) | 333 s | 14.1 GB |
| `count3` (compact, adapter-masked) | 95 organisms: 23 above, 50 from a run's organism list, 11 microbes, 16 class representatives, RefSeq viruses | ~119 Gb | **47.9 GB** (43.6 GB entries + 4.3 GB table) | 49 min (`--jobs 3`) | 13.6 GB |

- Index size grows **linearly** with total genome size: 0.06–0.12 minimizers per base.
  With 4-byte entries that is **0.35 bytes per base** plus the fixed 4.3 GB table
  (0.49 bytes per base overall for `all23c`); with 8-byte entries it was 0.71 bytes per base
  (2.77 G entries for 31.4 Gb). `count3` measured 0.37 bytes per base of entries.
- Repeat-rich genomes cost less per base, because a repeated k-mer is stored once:
  common toad 0.058 and nettle 0.064 minimizers per base, against about 0.10 for mammals.
- Shared entries are rare between different species: 35 M of 2.77 G entries (1.3 %).
  Every genome in `all23` keeps ≥ 90 % unique markers (the two songbirds are lowest, at 0.91).
  Not so within a family: in `count3`, grey and ringed seal kept 0.11, goat, sheep and muskox
  0.30–0.32, four birds-of-paradise 0.20–0.26; grouped under one name each they keep 0.66–0.88.
  Pairs such as cat and lion, dog and arctic fox, moose and reindeer still keep only 0.46–0.51.
- Peak build RAM depends on `--jobs` × the largest genomes being processed at once,
  not on the total. The build also writes temporary run files (`<out>.tmpN`) about as
  large as the final index, and deletes them at the end.

### Scan

One paired-end lane pair (lane A) from a NovaSeq run of 246 samples (the benchmark run) (~0.37 GB gz R1+R2, ~660 Mb),
index warm in page cache:

| Index | `--every 50` | `--every 1` | Scan RSS (mostly shared page cache) |
|---|---|---|---|
| `all23` (8-byte, binary search over the whole array) | 7.2 s | 196 s | 21.9 GB |
| `all23c` (compact, prefix table) | **3.5 s** | **53 s** | 15.2 GB |

`count3` (47.9 GB), R1 only of lane A (2.22 M reads), `-t 8`, index in
`/var/tmp`: **1.9 s** at `--every 50`, **37.6 s** at `--every 1`, 47.8 GB RSS. The same scans
with the index in the ecryptfs home took 43 s and 966 s (see Hard limits).

The 16 lane pairs used for benchmarking gave byte-identical reports with both indexes.
The earlier runs below used the 8-byte index:

| Workload | Time | Notes |
|---|---|---|
| 16 lane pairs, 8 at a time, `--every 50` | 18.6 s | 348 MB/s of gzipped input overall |
| Whole project folder: 492 lane pairs, 181 GB gz, 8 at a time | 504 s (8.4 min) | no errors |
| Same folder with `bycatch run -j 8 -t 2` (compact index, current code) | 256 s (4.3 min) | no errors; 17.9 GB peak RSS, mostly the shared index |
| Single uncompressed R1, `--every 50` | ~630 MB/s | about 12–15 min per 100 GB gzipped file |
| 100 Mb of simulated long reads, `--every 1` | 16 s | mostly index lookups |

- **Memory per scan**: 340 MB private. Almost all of it is the "seen" bitmap, one bit
  per index entry (2.77 G entries / 8). The index itself (15.4 GB compact) is page
  cache shared by all scans running at once. `top` therefore shows ~15 GB RSS per
  process, but 8 parallel scans still need only ~15 GB + 8 × 340 MB. With `count3`
  (10.9 G entries) the bitmap is 1.4 GB per scan: keep `run -j` at about 4 on this machine.
- The first scan after a reboot reads the index from disk once. Later scans find it
  in page cache.

### Accuracy on simulated data (10-genome index)

- Short reads (150 bp, 0.5 % errors) and long reads: no read assigned to the wrong organism.
- Spiked contamination: estimated 10.54 % against 10.55 % true for short reads; exact for long reads.
- `--every 50` on 500,000 reads gave 10.35 %, within the expected sampling noise (±0.3 % at 10,000 sampled reads).

Fraction of human reads classified, by read length (w = 15):

| Read length | Classified |
|---|---|
| 150 nt | 99.8 % |
| 100 nt | 99.3 % |
| 50 nt | 65 % |

## 2. Caveats (how to read the results)

1. **Percentages are of fragments, not of bases or cells.** For long reads a 100 kb read
   and a 1 kb read count the same. For paired reads one pair counts once.
2. **Organisms not in the index do not disappear.** Their reads show as `unmatched`, or as
   false hits on conserved sequence (16S/23S rRNA, mitochondria, ubiquitous genes) to the
   closest organism in the index. Check the `diversity` column: those hits pile onto few
   markers and get the `low_diversity` note (real genomes measured 0.59–1.0, 16S amplicons 0.005–0.16).
3. **`low_diversity` can also flag a real contaminant** with heavy PCR duplication or very few
   reads, because duplicates hit the same markers again. The note means "check this", not "artefact".
4. **Gaps in the host assembly look like contamination or `unmatched`.** For example, a trace
   of a second organism (0.03–0.04 %) found in every sample of one run may be real DNA
   of that organism or host sequence missing from the host assembly; the reads alone
   cannot tell the two apart. Very low signals (< 0.1 %)
   in every sample of a run should be read with this in mind. `bycatch run` flags such
   organisms as `background` (in ≥ 80 % of samples at a similar level).
5. **Close relatives cancel each other.** A k-mer present in two genomes becomes `shared` and
   counts for neither. Two strains of one species, or two species of one genus, lose most of
   their markers. `build` warns when an organism keeps less than 50 % unique markers.
   Group such genomes under one name instead.
6. **Exact 27-mer matching needs close sequence identity.** A contaminant that is a different
   strain or species from the reference (more than a few % divergent) gives fewer hits and may
   end up as `weak` or `unmatched`. High-error long reads (ONT, ~5–10 %) also lose k-mers, but
   their length usually makes up for it.
7. **Sampling noise.** With `--every N` the standard error of a fraction p is
   √(p(1−p)/n), where n is the number of sampled fragments. At n = 100,000 a 0.04 % contaminant
   is measured from ~40 reads (±16 %). Low-level results from small files need a smaller N.
   The report gives 95 % Wilson intervals (`ci95_low`, `ci95_high`, `# contamination_ci95`);
   they cover this sampling noise and duplicates (design effect D, about 1.00 at
   `--every 50` and 1.09 at `--every 1` on the benchmark run's lanes), not classification bias.
   With `--max-reads` or `--precision` the sampled reads come only from the start of the
   file (the first tiles), so the result may differ slightly from a sample spread over the
   whole file.
8. **Adapters and primer dimers** are reported as `adapter`, poly-G and simple repeats as
   `low_complexity` (exact 17-mers of a built-in Illumina/MGI list; other kits stay
   `unmatched`). On one lane of the benchmark run they were 0.55 % and 0.19 % of fragments; in the worst
   one 2.6 % and 0.4 %, next to 16 % `unmatched`. `--unmatched-out` writes a sample of the
   rest for BLAST or sylph. Reference genomes are adapter-masked at build; before that,
   adapter remnants in RefSeq phages and *Picea abies* turned 0.2 % of one lane into
   `viruses` hits and raised `conflict` from 0.8 % to 3.5 %. About 0.2 % of that lane's reads
   still carry adapter but hit a genome, likely adapter variants off the built-in list.
9. **Per-read output (`--reads`)** writes the R1 sequence only, even for paired input.

### Run-level flags on the benchmark run (492 lane pairs, host cod, `--every 50`)

| Organism | Flag | Lane pairs | Reading |
|---|---|---|---|
| mouse | `background` | 491 of 492 | ~0.02 % in every sample: the trace of caveat 4 |
| common frog | `background` | 434 | same pattern |
| *E. coli* | `elevated` | sample B (0.34–0.37 %), sample C (0.06–0.07 %), both lanes | run median 0 % |
| human | `elevated` | sample D, both lanes (0.04–0.05 %) | run median 0.002 % |

Single-lane `elevated` calls of 5–6 reads (reed warbler, honeybee, *M. tuberculosis*) are at
the edge of the rule; a call seen in both lanes of a sample is the stronger signal. No
`index_hopping` flag fired: no sample had ≥ 1 % of a non-host organism.

## 3. Hard limits

| Limit | Value | Cause |
|---|---|---|
| Organisms per index | 255 | 8-bit organism id; 255 is the `shared` tag |
| Index size | must fit in free RAM (page cache) | lookups are random; from disk the scan is orders of magnitude slower |
| On this machine | ~55 GB free → ~140 Gb of genomes | 4.3 GB table + 0.35 bytes per base; `count3` (47.9 GB) uses most of it |
| Encrypted home (ecryptfs) | index must live elsewhere | ecryptfs caches each page twice (encrypted and decrypted), so the index needs twice its size in RAM; `count3` scans were 25× slower from there |
| Read length | ≥ 41 bp for at least one minimizer, ≥ ~100 bp for reliable calls | w = 15, k = 27, 2-hit rule |
| k-mer length | fixed at 27 | `kmer << 8 \| id` must fit in 64 bits |
| Paired files | R1 and R2 must have the same number of records in the same order | checked; a mismatch is an error |
| Classification threads | `--threads`, default 4 | plus one reader thread per input file; decompression stays single-threaded per file |
| Prefix-table offsets | < 2³² entries per block of 2¹⁶ prefixes | never reached in practice; `build` checks it |

## 4. Ways to reduce index size, memory and runtime

### Index size and RAM

| Strategy | Effect on `all23c` (15.4 GB) | Cost | Recommendation |
|---|---|---|---|
| **Drop `shared` entries** from the index | −1.3 % | `shared` reads become `unmatched` | Not worth it |
| **Group strains / close relatives** under one name | keeps markers unique; no extra space | cannot tell grouped organisms apart | Use when needed |
| **Leave out organelles and plasmids** | small | lose mitochondrial/chloroplast detection | Optional; also reduces 16S-like false hits from chloroplasts |
| **Larger window, w = 31** | entries halve, table stays: ~15.4 → ~9.8 GB (−36 %) | 50 nt reads become unusable (~1.5 minimizers per read); 150 nt reads keep ~8 | Only for libraries ≥ 150 nt |
| **Repeat masking** (RepeatMasker lowercase → N), measured on human | 306 M → 158 M minimizers (−48 %) | sensitivity 99.98 % → 52.3 %, uneven between species | Do not use |
| **Low-complexity masking (sdust)** | small (not measured) | little; `build` already drops k-mers with < 8 distinct 3-mers | Not worth it |
| **FracMinHash sketch** (sylph / sourmash style, s = 200) | 15.4 GB → ~1.3 GB (s = 1000: ~250 MB) | no per-read output; abundance estimates only | Best for hundreds of species or 50 nt reads |
| **Two-stage screen**: FracMinHash over many genomes, then the current classifier on the few organisms found | small stage-2 index | two steps | Best long-term design if the reference list grows |
| Alignment (FastQ Screen, BBSplit) | larger (Bowtie2 ~1.2 bytes per base) | 10–100× slower | Only as an independent cross-check |

### Scan memory

- The per-scan 340 MB bitmap scales with the number of index entries (not bytes), so
  the compact layout does not change it. A sparse set of seen positions would use less when few reads are sampled,
  but 340 MB is not a problem today.
- Running many scans in parallel costs one copy of the index in page cache, not one per scan.

### Runtime

- **At `--every 50` decompression sets most of the speed** (with the compact index). Options, in order of effect:
  - run several files in parallel (8 at a time in the 492-pair run above);
  - decompress with a parallel tool and pipe in: `rapidgzip -dc R1.fastq.gz | bycatch scan ... -`
    (only for single-end; paired input needs two files);
  - the `zlib-rs` backend is already in use.
- **At `--every 1` index lookups set the speed** (53 s against 3.5 s per lane pair with the
  compact index and one thread; the prefix table already cut it from 196 s). Classification
  now runs on `--threads` workers. Measured with the 10-genome index, `--every 1`, warm cache:

  | Input | before (1 thread) | `-t 1` | `-t 4` | `-t 16` |
  |---|---|---|---|---|
  | 500,000 simulated 150 nt reads | 35.9 s | 16.6 s | 4.3 s | 2.2 s |
  | 20,000 simulated long reads | 21.9 s | 21.0 s | 5.3 s | 2.8 s |

  | One real lane pair (lane A, 2.22 M pairs, `all23c`) | 34.1 s | 30.4 s | 9.3 s | 5.1 s (`-t 8`: 6.0 s) |

  `-t 1` is faster than before because decompression now runs on its own thread. Above
  8 threads the gain flattens: R1 decompression alone takes 2.7 s, and random lookups
  compete for memory bandwidth. Reports and `--reads` files were byte-identical to the old
  ones (apart from the new columns) and between thread counts; so were the reports of all
  16 benchmark lane pairs at `--every 50`.
  `classify` now prefetches the table slots of all of a fragment's minimizers, then the
  first entry of each range, before searching any, so the cache misses overlap. With
  fat LTO (`[profile.release]`) on top, measured on `count3`, lane A, `-t 8`, warm
  cache, best of 3: R1 `--every 1` 4.2 → 3.3 s, the pair `--every 1` 8.3 → 6.4 s;
  `--every 50` unchanged at 1.0 s (decompression bound). Reports and `--reads` files
  byte-identical. Tried and dropped, both measured slower: a plain load before the
  bitmap's `fetch_or`, and a calloc'd (lazily zeroed) bitmap.
- **Build time** (5.4 min for 31 Gb) is dominated by parsing and sorting each genome.
  `--jobs` above 4 helps if RAM allows. Adding genomes does not need a rebuild:
  `bycatch add` reads only the new genomes plus one pass over the old index (decoding
  `all23c`'s 2.77 G entries takes about 70 s).
