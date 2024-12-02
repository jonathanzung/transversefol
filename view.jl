using Serialization

include("batch/2cusp_manifest.txt")
include("search.jl")

isosig="ivvPQQcfghghfhgfaddddaaaa_20000222"

tup = deserialize("/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls")
bt=tup.bt
Eupper=tup.Eupper
Elower=tup.Elower
Elong=tup.Elong

p=PlotlyJS.plot()
addtraces!(p, _plotjs(Elower, Eupper)...)
add_trace!(p, _plotjs(Elong;color=LONGITUDE_COLOUR))

display(p)
