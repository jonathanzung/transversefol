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

veering_census = CSV.read(joinpath(@__DIR__,"..","data","veering_census_with_data.txt"), DataFrame,
             types = [String, String, Int,  String, Int, String, Int, Vector{Int}, Vector{Int}, String, Vector{String}],
             header = [:isosig, :depth, :ncusps, :geometric, :nsymmetries, :edgeorientable, :eulerclass, :nladders, :ntetrahedra, :homology, :names]
            )

const by_cusps = [subset(veering_census, :ncusps => (n->n.==i)) for i in 1:6]

const depth1 = subset(veering_census, :ncusps => (n->n.==2), :depth => (n-> n.=="F1"))

const depth0 = subset(veering_census, :ncusps => (n->n.==2), :depth => (n-> n.=="F0"))

function index(isosig::String)
    for i in 1:nrow(veering_census)
        if veering_census[i, :isosig] == isosig
            return i
        end
    end
    return 0
end
function lookup(i)
    return veering_census[i,:isosig]
end
function lookup(i, ncusps)
    return by_cusps[ncusps][i,:isosig]
end

function lookup_row(i)
    return veering_census[i,:]
end
function lookuprow(i, ncusps)
    return by_cusps[ncusps][i,:]
end

end
