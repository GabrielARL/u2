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
yb = signal(repeat(mseq(12); inner=12), 6000)
filename = "loc6-1.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
bbc , bb2 = prepsig(x[:,2])

 #events = Int[]
 #events = find_mseq(bbc, (yb), 0.8, 3000) #0.7, 3000 for 4-2 #0.5, 3000 for 5-1 # 0.7, 2000 for 5-2 # 0.7,2000
 #pushfirst!(events, 228395)
 #pushfirst!(events, 98365)

#timeoffset_4_2 = [0 23 0 920 687 2158 717 263]'
#events = events .+  timeoffset_4_2


#events = [105659, 224954, 371316, 471439, 563463+30000-989, 718312-177, 845395-4000+83, 966890+31000]

cber1s = Float64[]
cber2s = Float64[] 
cber3s = Float64[] 
cber4s = Float64[] 

for i = 1:length(events)
 cber1, cber2, cber3 , cber_com = process3(x, events[i], 10, 10)
 push!(cber1s, cber1)
 push!(cber2s, cber2)
 push!(cber3s, cber3)
 push!(cber4s, cber_com)
 println((cber1,cber2,cber3,cber_com))
end

save("loc6.jld2", "cber1s", cber1s, "cber2s", cber2s, "cber3s", cber3s, "cber4s", cber4s)