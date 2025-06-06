using Base.Threads
using VeeringCensus

include("batch/1cusp_manifest.txt")

isosigs=[row.isosig for row in eachrow(VeeringCensus.depth1)][1:10]
@threads for l in isosigs
    println(l)
    run(`python3 prepare.py $(l)`)
end
