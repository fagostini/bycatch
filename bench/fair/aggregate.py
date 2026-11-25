#!/usr/bin/env python3
"""Aggregate count tables into result tables, using the frozen score.py and the pre-registered
threshold rule. Usage: aggregate.py counts_dir out_dir"""
import csv, glob, os, re, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score import COLS, score  # noqa: E402

CNT, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
SETRE = re.compile(r"(?:(?P<setting>[a-z]+[0-9.]+)_)?L(?P<L>\d+)_(?P<mix>low|mid|high)_s(?P<seed>\d)$")
DEFAULT = {"kraken2_ncbi": "c0", "kraken2_flat": "c0", "kunpeng_ncbi": "c0", "kunpeng_flat": "c0",
           "deacon": "a2", "bycatch": "-", "sourmash": "th50000"}


def truth(set_name):
    rows = [l.rstrip("\n").split("\t") for l in open(f"{CNT}/truth/{set_name}.truth")]
    return [(o, int(k), c) for o, k, c in rows]


rows = []
for tool in sorted(d for d in os.listdir(CNT) if d not in ("truth", "sylph", "sourmash")):
    for f in sorted(glob.glob(f"{CNT}/{tool}/*.tsv")):
        m = SETRE.match(os.path.basename(f)[:-4])
        counts = {}
        for line in open(f):
            t, c, v = line.rstrip("\n").split("\t")
            counts[(t, c)] = int(v)
        set_name = f"L{m['L']}_{m['mix']}_s{m['seed']}"
        r = score(counts, truth(set_name))
        rows.append({"tool": tool, "setting": m["setting"] or "-", "L": int(m["L"]), "mix": m["mix"],
                     "seed": int(m["seed"]), **r})

# sample-level tools
for f in sorted(glob.glob(f"{CNT}/sylph/profile_*.tsv")) + sorted(glob.glob(f"{CNT}/sourmash/*.csv")):
    base = os.path.basename(f)
    tool = "sylph" if base.startswith("profile_") else "sourmash"
    setting = "profile" if tool == "sylph" else base.split("_")[0]
    set_name = base.split("_", 1)[1].rsplit(".", 1)[0]
    m = SETRE.match(set_name)
    tr = truth(set_name)
    tot = sum(k for _, k, _ in tr)
    inp = sum(k for _, k, c in tr if c == "in_panel") / tot
    est, found = 0.0, set()
    with open(f) as fh:
        rd = csv.DictReader(fh, delimiter="\t" if tool == "sylph" else ",")
        n = 0
        for x in rd:
            if tool == "sylph" and not x.get("Genome_file"):
                continue  # timing line appended by T, not a sylph row
            n += 1
            if tool == "sylph":
                name = os.path.basename(x["Genome_file"])
                frac = float(x["Sequence_abundance"]) / 100
            else:
                name, frac = x["name"], float(x["f_unique_weighted"])
            org = name
            found.add(org)
            if "cod" not in org.lower():
                est += frac
    if tool == "sylph" and n == 0:
        continue  # profile hit the 1800 s timeout: no output, reported as missing
    rows.append({"tool": tool, "setting": setting, "L": int(m["L"]), "mix": m["mix"], "seed": int(m["seed"]),
                 "inpanel_est": est, "inpanel_true": inp, "inpanel_abs_err": abs(est - inp),
                 "detected": len(found) - any("cod" in x.lower() for x in found)})

keys = ["tool", "setting", "L", "mix", "seed"] + COLS
with open(f"{OUT}/all_runs.tsv", "w") as o:
    o.write("\t".join(keys) + "\n")
    for r in rows:
        o.write("\t".join(str(r.get(k, "")) if not isinstance(r.get(k), float) else f"{r[k]:.6g}" for k in keys) + "\n")

# pre-registered selection on L150 mid seed 1
sel = {}
for tool in {r["tool"] for r in rows}:
    cand = [r for r in rows if r["tool"] == tool and r["L"] == 150 and r["mix"] == "mid" and r["seed"] == 1]
    if not cand:
        continue
    if tool.startswith(("kraken2", "kunpeng")):
        key, best = "macro_f1", max
    elif tool == "deacon":
        key, best = "contam_f1", max
    else:
        key, best = "inpanel_abs_err", min
    target = best(r[key] for r in cand)
    winners = [r["setting"] for r in cand if r[key] == target]
    sel[tool] = DEFAULT.get(tool) if DEFAULT.get(tool) in winners else winners[0]
with open(f"{OUT}/selected.tsv", "w") as o:
    for t, s in sorted(sel.items()):
        o.write(f"{t}\t{s}\t(default {DEFAULT.get(t, '-')})\n")

# reported numbers: seeds 2-3, mean, at the selected and default settings
with open(f"{OUT}/summary.tsv", "w") as o:
    o.write("\t".join(["tool", "setting", "which", "L", "mix"] + COLS) + "\n")
    for tool in sorted(sel):
        for which, s in (("selected", sel[tool]), ("default", DEFAULT.get(tool, sel[tool]))):
            if which == "default" and s == sel[tool]:
                continue
            for L in (50, 100, 150):
                for mix in ("low", "mid", "high"):
                    rs = [r for r in rows if r["tool"] == tool and r["setting"] == s and r["L"] == L
                          and r["mix"] == mix and r["seed"] in (2, 3)]
                    if not rs:
                        continue
                    vals = []
                    for c in COLS:
                        xs = [r[c] for r in rs if c in r and r[c] == r[c]]
                        vals.append(f"{statistics.mean(xs):.6g}" if xs else "")
                    o.write("\t".join([tool, s, which, str(L), mix] + vals) + "\n")
print(open(f"{OUT}/selected.tsv").read())
