module VeeringCensus

using CSV
using DataFrames
import Base: tryparse
export veering_census

function tryparse(::Type{Vector{Int}}, str::String)
    parse.(Int, split(chop(str, head=1),','))
end
function tryparse(::Type{Vector{String}}, str::String)
    split(chop(str, head=1),',')
end

veering_census = CSV.read("veering_census_with_data.txt", DataFrame,
             types = [String, String, Int,  String, Int, String, Int, Vector{Int}, Vector{Int}, String, Vector{String}],
             header = [:isosig, :depth, :ncusps, :geometric, :nsymmetries, :edgeorientable, :eulerclass, :nladders, :ntetrahedra, :homology, :names]
            )

onecusp = subset(veering_census, :ncusps => (n->n.==1))
twocusp = subset(veering_census, :ncusps => (n->n.==2))
threecusp = subset(veering_census, :ncusps => (n->n.==3))

depth1 = subset(veering_census, :ncusps => (n->n.==2), :depth => (n-> n.=="F1"))

end
