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
wavsize(filename) = wavread(filename; format="size")

y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
filename = "loc5-1.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
bbc, bb2 = prepsig(x[:,2])

# events = Int[]
# events = find_mseq(bbc, (yb), 0.7, 2000) #0.7, 3000 for 4-2 #0.5, 3000 for 5-1 # 0.7, 2000 for 5-2 # 0.7,2000
# e_all = Float64[]


snr1,snr2,snr3 = process3(x, events[2], 1, 4)

#uber, cber = process(bbc, events[1], 1, 4);