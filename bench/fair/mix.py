#!/usr/bin/env python3
"""Mix ART pools into one synthetic paired set with exact, known proportions.
Usage: mix.py L seed mix outprefix      (mix: low | mid | high | thru)
Writes outprefix_R1.fq.gz, outprefix_R2.fq.gz, outprefix.truth (organism, pairs, class).
Read names are <organism>_<n>: the truth for every read."""
import random, subprocess, sys

POOLS = "/var/tmp/bench_fair/pools"
IN_PANEL = ["human", "mouse", "lsinapis", "drosophila", "arabidopsis", "yeast", "mtb", "ecoli", "phix"]
OFF_PANEL = ["pacific_cod", "dog", "frog"]
MIXES = {  # pairs, in-panel fraction, off-panel fraction; the rest is cod
    "low": (2_000_000, 0.0007, 0.0003),
    "mid": (500_000, 0.007, 0.003),
    "high": (500_000, 0.07, 0.03),
    "thru": (10_000_000, 0.007, 0.003),
}
AD1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # TruSeq read 1 adapter
AD2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # TruSeq read 2 adapter
POLYG_RATE = 0.02  # NovaSeq dark-cycle tails


def split(total, names):
    base, rem = divmod(total, len(names))
    return {n: base + (i < rem) for i, n in enumerate(names)}


def read_pool(path):
    with open(path) as f:
        lines = f.read().split("\n")
    return [(lines[i + 1], lines[i + 3]) for i in range(0, len(lines) - 3, 4)]


def damage(seq, qual, L, adapter, insert, rng):
    if insert < L:  # adapter read-through: the insert ends inside the read
        tail = (adapter + "G" * L)[: L - insert]
        seq, qual = seq[:insert] + tail, qual[:insert] + qual[insert:][: len(tail)]
    if rng.random() < POLYG_RATE:
        start = rng.randrange(L // 2, L)
        seq, qual = seq[:start] + "G" * (L - start), qual[:start] + "#" * (L - start)
    return seq, qual


def main():
    L, seed, mix, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    n, fin, foff = MIXES[mix]
    pool_seed = 9 if mix == "thru" else seed
    counts = split(round(n * fin), IN_PANEL) | split(round(n * foff), OFF_PANEL)
    counts["cod"] = n - sum(counts.values())
    rng = random.Random(f"{L}-{seed}-{mix}")
    pairs = []
    for org, k in counts.items():
        r1 = read_pool(f"{POOLS}/{org}_L{L}_s{pool_seed}_1.fq")
        r2 = read_pool(f"{POOLS}/{org}_L{L}_s{pool_seed}_2.fq")
        assert len(r1) == len(r2) >= k, f"{org}: pool {len(r1)} < {k}"
        for j in rng.sample(range(len(r1)), k):
            pairs.append((org, r1[j], r2[j]))
        del r1, r2
    rng.shuffle(pairs)
    # shortcut: R1/R2 tails are cut independently of the true fragment ends; fine for k-mer tools.
    g1 = subprocess.Popen(["gzip", "-1"], stdin=subprocess.PIPE, stdout=open(f"{out}_R1.fq.gz", "wb"), text=True)
    g2 = subprocess.Popen(["gzip", "-1"], stdin=subprocess.PIPE, stdout=open(f"{out}_R2.fq.gz", "wb"), text=True)
    for i, (org, (s1, q1), (s2, q2)) in enumerate(pairs):
        ins = max(20, round(rng.gauss(250, 80)))
        s1, q1 = damage(s1, q1, L, AD1, ins, rng)
        s2, q2 = damage(s2, q2, L, AD2, ins, rng)
        name = f"{org}_{i}"
        g1.stdin.write(f"@{name}\n{s1}\n+\n{q1}\n")
        g2.stdin.write(f"@{name}\n{s2}\n+\n{q2}\n")
    for g in (g1, g2):
        g.stdin.close()
        assert g.wait() == 0
    with open(f"{out}.truth", "w") as t:
        for org, k in counts.items():
            cls = "host" if org == "cod" else ("in_panel" if org in IN_PANEL else "off_panel")
            t.write(f"{org}\t{k}\t{cls}\n")


if __name__ == "__main__":
    main()
