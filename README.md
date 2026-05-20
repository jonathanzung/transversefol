# Dependencies
Requires a python3 with snappy, regina, and veering (currently requires master branch on github)

# Installation
Set your python path with the `TRANSVERSEFOL_PYTHON` environment variable.

Make two directories (used by TransverseFol for cache), and point to them with the environment variables `TRANSVERSEFOL_CACHE_DIR` and `TRANSVERSEFOL_PREP_CACHE_DIR`.

In the Julia terminal, run

```
julia> using Pkg
julia> Pkg.develop(path="/path/to/repo")
julia> Pkg.resolve()
```

# Usage

Start julia with your desired number of threads (8):
```
bash> julia -t 8 
```

Load transversefol and a veering triangulation
```
julia> using Revise #recommended; install if needed
julia> using TransverseFol
julia> phi = TransverseFol.load("eLMkbcddddedde_2100")
```

Search for pseudo-Anosov flows given a manifold name. Works with cusped manifolds as well. `scripts` is a top level directory in the repo.
```
python3 scripts/find_pA_flows.py "m304(4,-1)"
```

Find foliations transverse to phi. Results will be cached in `TRANSVERSEFOL_CACHE_DIR`. If you didn't find enough foliations, you can try to run it again, or replace `TRY` with `TRYHARD`
```
julia> runjob(phi; TRY...) 
julia> quickview(phi)
```

Search for foliations on a closed 3-manifold by searching for pseudo-Anosov flows, and then searching for transverse foliations
```
julia> include("scripts/find_foliation.jl") 
julia> find_foliation("m304(4,1)"; TRY...)
```

# Ziggurat zoo 
See some sample output here: https://web.mit.edu/jzung/www/ziggurat_zoo/ziggurat_zoo.html
