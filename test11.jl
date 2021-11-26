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

filename = "loc4-2.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
bbc=prepsig(x[:,2])
events = Int[]
events = find_mseq(bbc, (yb))
T = length(yb)
rbb = bbc[events[1]:events[1]+T]
vn = makevb(rbb) 
vn1 = vn#resample(vn[1:min(length(yb), length(vn))],(1000)/1000)
blksize = 1200

e = Float64[]
e2 = Float64[]
SDs = ComplexF64[]
atxDs = ComplexF64[]

vn_best = zeros(ComplexF64,1,blksize)
aSD_zm = zeros(ComplexF64,1,blksize)
tx_zm = zeros(ComplexF64,1,blksize)

for i = 1:blksize:length(vn1)-blksize
    vn_best = compDoppler(vn1[i:i+blksize], yb[i:i+blksize])
    if(length(vn_best) > blksize)
        vn_best = vn_best[1:blksize]
    end
    if(length(vn_best) < blksize)
        vn_best = [vn_best; zeros(ComplexF64,blksize-length(vn_best))]
    end
    quantizer_op, error, ff_filt, SD, ffops, ff = chRlsTrain(vn_best[1:blksize], yb[i:i+blksize], 4*12, 0.9995, 1)
    quantizer_op, SD2 = DA(vn_best[1:blksize], ff_filt, 1)
    atx=angle.(samples(yb[i:i+blksize]))
    tx_m=mean(atx)
    aSD = abs.(angle.(SD2))
    tx_zm = atx .- tx_m
    aSD_zm = aSD .- tx_m
    println(length(aSD_zm))
    downbb, downtxbb = downsymboling(aSD_zm, tx_zm)
    ber = computeBER(downbb, downtxbb)
    println((i,ber))
    push!(e2,ber)
end







