using Revise
using JSON
includet("script.jl")

#@show JSON.parse("batch/$(VeeringCensus.lookup(1,2)).json")


tup=load("nLLvLMAMQkbcijikhijlkmmmtsnjnasfnihbkw_2010222112011",refresh_prep=true, nlongs=1)
cand=tup.Elong.A[1][2]
display(draw(cand,1,curve=tup.bt.snappy_weights))
