#!/usr/bin/env python3
"""Turn one tool's per-read output into a count table (true_organism, call, pairs). Frozen before runs.

Usage:
  calls.py bycatch  reads.tsv                       # bycatch scan --reads output
  calls.py kraken   out.txt  taxmap  truth.tsv      # Kraken2 / Kun-peng per-read output
  calls.py deacon   host_names contam_names truth.tsv
The true organism is the read-name prefix (<organism>_<n>). Pairs missing from an output
(Kun-peng omits unclassified reads) are counted as "unassigned" via the truth totals.
taxmap: TSV  taxid  call  -- precomputed per DB: every taxid in the DB's taxonomy mapped to the panel
organism whose subtree holds it, else "higher" (rank above all panel species, holding >1) or
"unassigned" (root and taxa outside every panel subtree).
"""
import collections, sys

SPECIAL = {"shared", "weak", "conflict", "adapter", "low_complexity", "unmatched"}


def org(name):
    return name.rsplit("_", 1)[0]


def write(c):
    for (t, k), v in sorted(c.items()):
        print(f"{t}\t{k}\t{v}")


def fill_missing(c, truth):
    seen = collections.Counter()
    for (t, _), v in c.items():
        seen[t] += v
    for o, k, _ in truth:
        if k > seen[o]:
            c[(o, "unassigned")] += k - seen[o]


def truth_of(path):
    return [(o, int(k), cl) for o, k, cl in (l.rstrip("\n").split("\t") for l in open(path))]


def main():
    mode = sys.argv[1]
    c = collections.Counter()
    if mode == "bycatch":
        for line in open(sys.argv[2]):
            rid, cat = line.split("\t", 2)[:2]
            if rid in ("read_id", "read"):  # header
                continue
            c[(org(rid), "unassigned" if cat in SPECIAL else cat)] += 1
        write(c)
    elif mode == "kraken":
        taxmap = dict(l.rstrip("\n").split("\t") for l in open(sys.argv[3]))
        for line in open(sys.argv[2]):
            f = line.split("\t", 3)
            call = "unassigned" if f[0] == "U" else taxmap.get(f[2].split()[-1].rstrip(")"), "unassigned")
            c[(org(f[1]), call)] += 1
        fill_missing(c, truth_of(sys.argv[4]))
        write(c)
    elif mode == "deacon":
        host = {l.strip() for l in open(sys.argv[2])}
        contam = {l.strip() for l in open(sys.argv[3])}
        for r in host | contam:
            c[(org(r), "both" if r in host and r in contam else ("cod" if r in host else "contaminant"))] += 1
        fill_missing(c, truth_of(sys.argv[4]))
        write(c)


if __name__ == "__main__":
    main()
