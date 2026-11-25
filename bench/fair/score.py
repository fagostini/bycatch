#!/usr/bin/env python3
"""Score the fair benchmark from count tables. Frozen before any tool run (see PROTOCOL.md).

Count table (one per tool x setting x set): TSV  true_organism  call  pairs
  call = a panel organism name, "contaminant" (deacon: hit only the contaminant index),
         "both" (deacon: hit both indexes), "higher" (taxonomic tools: a rank above species
         that holds no single panel organism), or "unassigned".
Truth file (from mix.py): organism  pairs  class(host | in_panel | off_panel)

Usage: score.py counts.tsv truth.tsv  -> one TSV row of metrics on stdout
       score.py --header             -> the header
"""
import math, sys

PANEL = ["cod", "human", "mouse", "lsinapis", "drosophila", "arabidopsis", "yeast", "mtb", "ecoli", "phix"]
SCORED = [o for o in PANEL if o not in ("cod", "phix")]  # phiX: best case (indexed sequence), footnote only
DETECT_MIN = 10  # pairs called an organism for it to count as detected

COLS = ["host_fp", "host_fp_lo", "host_fp_hi", "host_lost", "macro_f1", "inpanel_est", "inpanel_true",
        "inpanel_abs_err", "nonhost_est", "nonhost_true", "nonhost_abs_err", "off_cod", "off_other_org",
        "off_unassigned", "off_flagged", "detected", "phix_recall", "contam_f1"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def score(counts, truth):
    n_true = {o: k for o, k, _ in truth}
    cls = {o: c for o, _, c in truth}
    total = sum(n_true.values())
    get = lambda t, c: counts.get((t, c), 0)
    calls_of = lambda c: sum(v for (t, cc), v in counts.items() if cc == c)
    host = n_true["cod"]
    # a host pair is a false positive when it is called contamination: any non-cod panel organism,
    # or deacon's contaminant-only call
    fp = sum(v for (t, c), v in counts.items() if t == "cod" and (c in PANEL[1:] or c == "contaminant"))
    lost = host - get("cod", "cod") - get("cod", "both")
    lo, hi = wilson(fp, host)
    named = any(c in SCORED for (_, c) in counts)  # deacon and other host/non-host tools name no organism
    f1s, detected = [], 0
    for o in SCORED:
        tp, called, real = get(o, o), calls_of(o), n_true.get(o, 0)
        p = tp / called if called else 0.0
        r = tp / real if real else 0.0
        f1s.append(2 * p * r / (p + r) if p + r else 0.0)
        detected += called >= DETECT_MIN
    contam_calls = sum(calls_of(o) for o in PANEL[1:]) + calls_of("contaminant")
    inpanel_true = sum(k for o, k, c in truth if c == "in_panel") / total
    nonhost_true = 1 - host / total
    cod_calls = calls_of("cod") + calls_of("both")
    off = [o for o, _, c in truth if c == "off_panel"]
    n_off = sum(n_true[o] for o in off)
    off_cod = sum(get(o, "cod") + get(o, "both") for o in off)
    off_org = sum(get(o, c) for o in off for c in PANEL[1:] + ["contaminant"])
    # binary: positive = an in-panel contaminant pair (phiX excluded); predicted = called contamination
    pos = [o for o in SCORED]
    is_contam = lambda c: c in PANEL[1:] or c == "contaminant"
    tp_b = sum(v for (t, c), v in counts.items() if t in pos and is_contam(c))
    pred_b = sum(v for (t, c), v in counts.items() if t != "phix" and is_contam(c))
    real_b = sum(n_true.get(o, 0) for o in pos)
    pb, rb = (tp_b / pred_b if pred_b else 0.0), (tp_b / real_b if real_b else 0.0)
    row = {
        "contam_f1": 2 * pb * rb / (pb + rb) if pb + rb else 0.0,
        "host_fp": fp / host, "host_fp_lo": lo, "host_fp_hi": hi, "host_lost": lost / host,
        "macro_f1": sum(f1s) / len(f1s) if named else float("nan"),
        "inpanel_est": contam_calls / total, "inpanel_true": inpanel_true,
        "inpanel_abs_err": abs(contam_calls / total - inpanel_true),
        "nonhost_est": 1 - cod_calls / total, "nonhost_true": nonhost_true,
        "nonhost_abs_err": abs(1 - cod_calls / total - nonhost_true),
        "off_cod": off_cod / n_off, "off_other_org": off_org / n_off,
        "off_unassigned": 1 - (off_cod + off_org) / n_off, "off_flagged": 1 - off_cod / n_off,
        "detected": detected if named else float("nan"),
        "phix_recall": get("phix", "phix") / n_true["phix"] if named else float("nan"),
    }
    return row


def main():
    if sys.argv[1] == "--header":
        print("\t".join(COLS))
        return
    counts = {}
    for line in open(sys.argv[1]):
        t, c, v = line.rstrip("\n").split("\t")
        counts[(t, c)] = counts.get((t, c), 0) + int(v)
    truth = [(o, int(k), c) for o, k, c in (l.rstrip("\n").split("\t") for l in open(sys.argv[2]))]
    assert sum(counts.values()) == sum(k for _, k, _ in truth), "counts do not sum to the input pairs"
    row = score(counts, truth)
    print("\t".join(f"{row[c]:.6g}" for c in COLS))


def _selftest():
    truth = [("cod", 990, "host"), ("human", 5, "in_panel"), ("mouse", 3, "in_panel"), ("phix", 1, "in_panel"),
             ("dog", 1, "off_panel")]
    counts = {("cod", "cod"): 985, ("cod", "human"): 2, ("cod", "unassigned"): 3, ("human", "human"): 5,
              ("mouse", "mouse"): 2, ("mouse", "unassigned"): 1, ("phix", "phix"): 1, ("dog", "unassigned"): 1}
    r = score(counts, truth)
    assert abs(r["host_fp"] - 2 / 990) < 1e-12 and abs(r["host_lost"] - 5 / 990) < 1e-12
    assert abs(r["inpanel_est"] - 10 / 1000) < 1e-12 and r["off_flagged"] == 1.0 and r["phix_recall"] == 1.0
    assert abs(r["nonhost_est"] - 15 / 1000) < 1e-12


if __name__ == "__main__":
    _selftest()
    main()
