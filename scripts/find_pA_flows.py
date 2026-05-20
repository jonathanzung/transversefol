"""
find_pA_flows.py — enumerate pA flows for a single closed manifold.

Usage:
    python3 find_pA_flows.py <manifold_name> [--method METHOD] [--count N]
                             [--max_drill N] [--maxlength F]

method: "combinatorial" (default) or "geodesic"

Prints one closed isosig per line, of the form:
    <isosig>_[<filling_slopes>]
"""

import sys
import os
import argparse
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pysrc.enumerate_pA import pA_flows
import snappy


def isom_to_json(isosig, cusp_images_MtoMM, cusp_maps_MtoMM):
    """Convert raw isometry data to BasisChange JSON format.

    Output keys:
      isosig    : the full "<base>_[<fillings>]" string
      slice     : [[p,q] per M cusp] filling slope or [0,0] if unfilled
      perm      : [MM cusp index (0-based)] for each surviving (unfilled) M cusp, in order
      cusp_maps : [2x2 matrix] one per surviving cusp, same order as perm
    """
    n = len(cusp_images_MtoMM)
    slice_ = []
    perm = []
    cusp_maps = []
    for j in range(n):
        img = cusp_images_MtoMM[j]
        m   = cusp_maps_MtoMM[j]
        if img is None:
            # Extract filling slope from the isosig string
            import ast
            idx = isosig.rfind("_[")
            fillings = ast.literal_eval(isosig[idx + 1:])
            slice_.append(list(fillings[j]))
        else:
            slice_.append([0, 0])
            perm.append(int(img))
            cusp_maps.append([[int(m[r, c]) for c in range(2)] for r in range(2)])
    return {"isosig": isosig, "slice": slice_, "perm": perm, "cusp_maps": cusp_maps}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("manifold_name")
    parser.add_argument("--method",       default="combinatorial")
    parser.add_argument("--count",        type=int,   default=10)
    parser.add_argument("--max_drill",    type=int,   default=3)
    parser.add_argument("--maxlength",    type=float, default=3.0)
    parser.add_argument("--max_segments", type=int,   default=6)
    parser.add_argument("--max_tets",     type=int,   default=20)
    parser.add_argument("--isometry",     action="store_true", default=False,
                        help="also print isometry data (slice, perm, cusp_maps)")
    parser.add_argument("--prong_counts", action="store_true", default=False,
                        help="also print prong counts for each filled cusp")
    args = parser.parse_args()

    # Parse manifold names like "m003(1,2)(3,4)" which snappy.Manifold() can't handle directly.
    # Split into base name and filling tuples.
    import re
    name = args.manifold_name
    filling_strs = re.findall(r'\(([^)]+)\)', name)
    base_name = re.split(r'\(', name)[0]
    if len(filling_strs) > 1:
        M = snappy.Manifold(base_name)
        fillings = [tuple(int(x) for x in s.split(',')) for s in filling_strs]
        M.dehn_fill(fillings)
    else:
        M = snappy.Manifold(name)
    need_extras = args.isometry or args.prong_counts
    for result in pA_flows(M, count=args.count, max_drill=args.max_drill,
                            maxlength=args.maxlength, max_segments=args.max_segments,
                            max_tets=args.max_tets, method=args.method,
                            return_isom=args.isometry,
                            return_prong_counts=args.prong_counts):
        if need_extras:
            s, extras = result
            obj = isom_to_json(s, *extras["isom"]) if args.isometry else {"isosig": s}
            if args.prong_counts:
                obj["prong_counts"] = extras["prong_counts"]
            sys.stdout.write(json.dumps(obj) + "\n")
        else:
            sys.stdout.write(result + "\n")
        sys.stdout.flush()
