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

function process(bbc, events)
    y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
    yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
    T = length(yb)
    rbb = bbc[events:events+T]
    vn = makevb(rbb) 
    vn1 = vn #resample(vn[1:min(length(yb), length(vn))],(1000)/1000)
    blksize = 1200
    e2 = Float64[]
    vn_best = zeros(ComplexF64,1,blksize)
    aSD_zm = zeros(ComplexF64,1,blksize)
    tx_zm = zeros(ComplexF64,1,blksize)
    push!(e2 , 0.0)
    length(vn1)
    for i = 1:blksize:length(vn1)-blksize
        vn_best = compDoppler(vn1[i:i+blksize], yb[i:i+blksize])
        if(length(vn_best) > blksize)
            vn_best = vn_best[1:blksize]
        end
        if(length(vn_best) < blksize)
            vn_best = [vn_best; zeros(ComplexF64,blksize-length(vn_best))]
        end
        quantizer_op, error, ff_filt, SD, ffops, ff = chRlsTrain(vn_best[1:blksize], yb[i:i+blksize], 4*12, 0.9995, 1)
        Ber1 = computeBER(quantizer_op, yb[i:i+blksize])
        #println((i,Ber1))
        quantizer_op, SD2 = DA(vn_best[1:blksize], ff_filt, 1)
        atx=angle.(samples(yb[i:i+blksize]))
        tx_m=mean(atx)
        aSD = abs.(angle.(SD2))
        tx_zm = atx .- tx_m
        aSD_zm = aSD .- tx_m
        #println(length(aSD_zm))
        downbb, downtxbb = downsymboling(aSD_zm, tx_zm)
        ber = computeBER(downbb, downtxbb)
        #println((i,ber))
        push!(e2,ber)
    end
    e2
end

y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
filename = "loc4-2.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
bbc, bbc2 = prepsig(x[:,2])
events = Int[]
events = find_mseq(bbc, (yb))
e_all = Float64[]
#events = events[2]
aaa = Float64[]
bbb = Int64[]
ccc = Int64[]

for i = 2:length(events)
    for j = 0:5:20 
        e2 = Float64[]
        temp = process(bbc, events[i]+j)
        println(length(temp))
        if(isassigned(temp))
            for v = 1:length(temp)
                push!(e2, temp[v])
            end
                if(!isempty(e2))
                    e2m = mean(e2)
                    push!(aaa,e2m)
                    push!(bbb,j)
                    push!(ccc,i)
                end
        end
    end    
        # for k = 1:length(temp)
        #     push!(e_all, temp[k])
        # end
end
#plot(e2)
#plot!(repeat([mean(e2)],length(e2)))





