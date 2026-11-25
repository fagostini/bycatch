#!/usr/bin/env python3
"""Map every taxid of a DB's taxonomy to a call for calls.py. Frozen before runs.
Usage: taxmap.py nodes.dmp panel_taxids.tsv > taxmap.tsv
panel_taxids.tsv: organism  taxid
Call: the panel organism whose subtree holds the taxid; "higher" for ancestors of a panel
organism (genus Gadus, Eukaryota, root, ...); "unassigned" for anything else."""
import sys

parent = {}
for line in open(sys.argv[1]):
    f = line.split("\t|\t")
    parent[f[0].strip()] = f[1].strip()
panel = dict(reversed(l.split()) for l in open(sys.argv[2]))  # taxid -> organism


def lineage(t):
    seen = [t]
    while parent.get(t, t) != t:
        t = parent[t]
        seen.append(t)
    return seen


ancestors = {a for t in panel for a in lineage(t)[1:]}
for t in parent:
    call = next((panel[a] for a in lineage(t) if a in panel), None)
    print(f"{t}\t{call or ('higher' if t in ancestors else 'unassigned')}")
