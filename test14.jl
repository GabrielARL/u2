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
using JLD2

include("rls.jl")
include("utils.jl")
include("FEC.jl")
include("BCH.jl")

import Main.FEC
import Main.BCH

wavsize(filename) = wavread(filename; format="size")

y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
filename = "loc5-1.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
bbc , bb2 = prepsig(x[:,2])

 #events = Int[]
 #events = find_mseq(bb2, (yb), 0.12, 5000) #0.7, 3000 for 4-2 #0.5, 3000 for 5-1 # 0.7, 2000 for 5-2 # 0.7,2000
 #e_all = Float64[]

events = [119417, 247465, 371316, 490000, 563463, 718312, 845395, 966890].+2000#events 5-1

 cber1s = Float64[]
 cber2s = Float64[] 
 cber3s = Float64[] 
 cber4s = Float64[] 

for i = 1:length(events)
 cber1, cber2, cber3 , cber_com = process3(x, events[i], 1)
 push!(cber1s, cber1)
 push!(cber2s, cber2)
 push!(cber3s, cber3)
 push!(cber4s, cber_com)
end

save("loc5.jld2", "cber1s", cber1s, "cber2s", cber2s, "cber3s", cber3s, "cber4s", cber4s)