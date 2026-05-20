"""
isosig_isometries.py — Given a list of veering isosigs with filling information,
find isometries to a common manifold (the fully-filled version of the first isosig)
and print JSON in the same format as find_pA_flows.py --isometry.

Usage:
    python3 isosig_isometries.py "base_[(p,q),(0,0)]" "base2_[(r,s),(0,0)]" ...

    or read from stdin (one isosig per line):
    cat isosigs.txt | python3 isosig_isometries.py

Each isosig has the form "<veering_isosig>_[<filling_slopes>]", e.g.:
    "eLMkbcddddedde_2100_[(2,-3),(0,0)]"

Cusps with slope (0,0) are left unfilled and correspond to cusps of the common
manifold MM.  MM is obtained by filling all non-(0,0) cusps of the first isosig.

Output (one JSON line per input isosig, same format as find_pA_flows --isometry):
    {"isosig": "...", "slice": [...], "perm": [...], "cusp_maps": [...]}

  slice[j]:     [p,q] filling slope for M cusp j, or [0,0] if unfilled
  perm[k]:      MM cusp index (0-based) for the k-th surviving M cusp
  cusp_maps[k]: 2x2 integer matrix M-SnaPPy -> MM-SnaPPy for the k-th surviving cusp
"""

import sys
import os
import ast
import json
import argparse
import random

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "find_pA"))
sys.path.insert(0, "/home/jonathan/Dropbox/repo/Veering/scripts")
sys.path.insert(0, "/home/jonathan/Dropbox/repo/Veering")

import snappy
import veering.taut


def parse_isosig(s):
    """Parse 'base_[fillings]' -> (base_isosig, fillings_list_of_tuples)."""
    s = s.strip()
    idx = s.rfind("_[")
    if idx == -1:
        return s, None
    base = s[:idx]
    fillings = ast.literal_eval(s[idx + 1:])
    return base, [tuple(f) for f in fillings]


def manifold_from_isosig(base_isosig):
    """Build a SnaPPy Manifold from a veering isosig, using veering orientation."""

    idx = base_isosig.rfind("_")

    if idx != -1: #then we're dealing with a veering isosig
        tri, _angle = veering.taut.isosig_to_tri_angle(base_isosig)
        assert tri.isOriented()
        return snappy.Manifold(tri)
    else: #it's a raw snappy isosig
        return snappy.Manifold(base_isosig)


def filled_manifold(M):
    """Return (MM, unfilled_indices) where MM has all Dehn-filled cusps removed."""
    is_complete = M.cusp_info('is_complete')
    filled_idx   = [i for i in range(M.num_cusps()) if not is_complete[i]]
    unfilled_idx = [i for i in range(M.num_cusps()) if     is_complete[i]]
    if not filled_idx:
        return M.copy(), unfilled_idx
    return snappy.Manifold(M.filled_triangulation(filled_idx)), unfilled_idx


def process_isosig(full_isosig, MM):
    """
    Find the isometry from M (with fillings applied, then filled) to MM.
    Returns a dict ready for JSON output, or None on failure.
    """
    base_isosig, fillings = parse_isosig(full_isosig)
    M = manifold_from_isosig(base_isosig)

    if fillings is None:
        fillings = [(0, 0)] * M.num_cusps()

    M.dehn_fill(fillings)
    M_filled, unfilled_idx = filled_manifold(M)

    try:
        isoms = M_filled.is_isometric_to(MM, return_isometries=True)
    except RuntimeError as e:
        print(f"# RuntimeError for {full_isosig}: {e}", file=sys.stderr)
        return None

    if not isoms:
        print(f"# No isometry found for {full_isosig}", file=sys.stderr)
        return None

    for isom in isoms:
        imgs = list(isom.cusp_images())
        maps = list(isom.cusp_maps())

        # Build BasisChange-style output: slice, perm, cusp_maps (no nulls).
        # unfilled_idx lists the surviving M cusps in order; for each, perm[k] = MM cusp.
        n = M.num_cusps()
        slice_ = []
        perm = []
        cusp_maps_out = []
        surviving_k = 0
        for j in range(n):
            if j in unfilled_idx:
                slice_.append([0, 0])
                perm.append(int(imgs[surviving_k]))
                m = maps[surviving_k]
                cusp_maps_out.append([[int(m[r, c]) for c in range(2)] for r in range(2)])
                surviving_k += 1
            else:
                f = fillings[j] if fillings else [0, 0]
                slice_.append(list(f))

        yield {
            "isosig":    full_isosig,
            "slice":     slice_,
            "perm":      perm,
            "cusp_maps": cusp_maps_out,
        }


def main():
    parser = argparse.ArgumentParser(
        description="Map veering isosigs with fillings to a common manifold and output isometry JSON.")
    parser.add_argument("isosigs", nargs="*",
                        help="Isosigs with filling info, e.g. 'base_[(p,q),(0,0)]'")
    args = parser.parse_args()

    isosigs = args.isosigs or [line.strip() for line in sys.stdin if line.strip()]
    if not isosigs:
        print("No isosigs provided.", file=sys.stderr)
        sys.exit(1)

    # Build MM from the first isosig
    first_base, first_fillings = parse_isosig(isosigs[0])
    M0 = manifold_from_isosig(first_base)
    if first_fillings is not None:
        M0.dehn_fill(first_fillings)
    MM, _ = filled_manifold(M0)
    print(f"# Common manifold MM: {MM}  ({MM.num_cusps()} cusp(s))", file=sys.stderr)

    for full_isosig in isosigs[1:]:
        for result in process_isosig(full_isosig, MM):
            sys.stdout.write(json.dumps(result) + "\n")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
