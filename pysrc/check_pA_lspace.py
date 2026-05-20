"""
check_pA_lspace.py — For each manifold in snappy.OrientableClosedCensus, check:
  - whether a pA flow was found (hodgson_weeks_pA/<name>_pAflows.txt is non-empty)
  - whether it is an L-space (from QHSpheres.csv, L_space column: 1=L-space, -1=not)

Falls back to matching via the 'descriptions' column if primary name not found.
Prints a summary table to stdout.
"""

import ast
import csv
import os
import pathlib

import snappy

CSV_PATH = "/home/jonathan/Downloads/conjecture_data/floer/final_data/QHSpheres.csv"
PA_DIR   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hodgson_weeks_pA")


def normalize(name):
    """Remove spaces: 'm003(-3, 1)' -> 'm003(-3,1)'."""
    return name.replace(" ", "")


# Build lookup from normalized name -> L_space value,
# including all aliases listed in the 'descriptions' column.
lspace_lookup = {}
with open(CSV_PATH, newline="") as f:
    for row in csv.DictReader(f):
        ls = row["L_space"]
        # primary name
        lspace_lookup[normalize(row["name"])] = ls
        # all aliases
        try:
            descs = ast.literal_eval(row["descriptions"])
            for d in descs:
                lspace_lookup[normalize(d)] = ls
        except Exception:
            pass


def has_pA_flow(normalized_name):
    """Returns True/False if computed, None if file absent."""
    fname = os.path.join(PA_DIR, normalized_name + "_pAflows.txt")
    p = pathlib.Path(fname)
    if not p.is_file():
        return None
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip()]
    return len(lines) > 0


print(f"{'name':<25} {'pA_flow':<10} {'L_space'}")
print("-" * 50)

counts = {"pA_lspace": 0, "pA_not_lspace": 0, "no_pA_lspace": 0, "no_pA_not_lspace": 0,
          "unknown_lspace": 0, "not_computed": 0}



L=[]

for M in snappy.OrientableClosedCensus(betti=0):
    name = str(M)
    key  = normalize(name)

    pA = has_pA_flow(key)
    ls = lspace_lookup.get(key)

    pA_str = "yes" if pA else ("no" if pA is False else "?")
    ls_str = "L-space" if ls == "1" else ("not L-space" if ls == "-1" else "unknown")

    print(f"{name:<25} {pA_str:<10} {ls_str}")

    if pA is None:
        counts["not_computed"] += 1
    elif ls is None:
        counts["unknown_lspace"] += 1
    elif pA and ls == "1":
        counts["pA_lspace"] += 1
    elif pA and ls == "-1":
        counts["pA_not_lspace"] += 1
    elif not pA and ls == "1":
        counts["no_pA_lspace"] += 1
    elif not pA and ls == "-1":
        L.append(key)
        counts["no_pA_not_lspace"] += 1

print()
print("Summary:")
for k, v in counts.items():
    print(f"  {k}: {v}")

print(L)