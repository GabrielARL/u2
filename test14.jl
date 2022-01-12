using Revise: push!, convert, copy
using ToeplitzMatrices: include
using LinearAlgebra: length
using Base: length_continued
using Statistics
using DataFrames
using Query
using LinearAlgebra
using CSV
using SignalAnalysis
using DSP
using WAV
using Plots
using Revise

include("rls.jl")
include("utils.jl")
include("FEC.jl")
include("BCH.jl")

import Main.FEC
import Main.BCH

wavsize(filename) = wavread(filename; format="size")

y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
filename = "loc6-1.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
bbc, bb2 = prepsig(x[:,2])

# events = Int[]
#  events = find_mseq(bbc, (yb), 0.7, 2000) #0.7, 3000 for 4-2 #0.5, 3000 for 5-1 # 0.7, 2000 for 5-2 # 0.7,2000
#  e_all = Float64[]

 aaa = Float64[]
 bbb = Float64[] 
 ccc = Float64[] 
 ddd = Float64[] 
 eee = Float64[] 
 fff = Float64[]
 ggg = Float64[]
for i = 1:1#length(events)
uber, cber = process(bbc, events[i], 1, 4);
 uber_com, cber_com, cber1, cber2, cber3 = process3(x, events[i], 1, 4)
 push!(aaa, uber_com)
 push!(bbb, cber_com)
 push!(ccc, cber1)
 push!(ddd, cber2)
 push!(eee, cber3)
 push!(fff, uber)
 push!(ggg, cber)
end

#save("loc6.jld2", "uber_com",aaa,"cber_com",bbb, "cber1", ccc, "cber2", ddd, "cber3", eee, "uber",fff,"cber",ggg)