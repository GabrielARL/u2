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

function process(bbc, events, d, Q)
    y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
    yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
    T = length(yb)
    rbb = bbc[events:events+T]
    vn = makevb(rbb, d) 
    vn1 = vn #resample(vn[1:min(length(yb), length(vn))],(1000)/1000)
    blksize = 300
    e2 = Float64[]
    vn_best = zeros(ComplexF64,1,blksize)
    aSD_zm = zeros(ComplexF64,1,blksize)
    tx_zm = zeros(ComplexF64,1,blksize)
    push!(e2 , 0.0)
    #length(vn1)
    λ = 0.9995 
    for i = 1:blksize:length(vn1)-blksize
        vn_best = compDoppler(vn1[i:i+blksize], yb[i:i+blksize])
        if(length(vn_best) > blksize)
            vn_best = vn_best[1:blksize]
        end
        if(length(vn_best) < blksize)
            vn_best = [vn_best; zeros(ComplexF64,blksize-length(vn_best))]
        end
        quantizer_op, error, ff_filt, SD, ffops, ff = chRlsTrain(vn_best[1:blksize], yb[i:i+blksize], Q*12, λ, 1)
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
        println((i,ber))
        push!(e2,ber)
    end
    e2
end





function process3(events, d, Q)
    y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
    yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
    T = length(yb)
    filename = "loc5-1.wav" 
    nsamples, _ = wavsize(filename)
    x = signal(filename; start = 1, nsamples = nsamples)
    bb1,_ = prepsig(x[:,1])
    bb2,_ = prepsig(x[:,2])
    bb3,_ = prepsig(x[:,3])
    println(events, length(bb1))
    rbb1 = bb1[events:events+T]
    rbb2 = bb2[events:events+T]
    rbb3 = bb3[events:events+T]

    vn1 = makevb(rbb1, d) 
    vn2 = makevb(rbb2, d)
    vn3 = makevb(rbb3, d)

    blksize = 300
    e2 = Float64[]
    vn_best = zeros(ComplexF64,1,blksize)
    aSD_zm = zeros(ComplexF64,1,blksize)
    tx_zm = zeros(ComplexF64,1,blksize)
    push!(e2 , 0.0)
    #length(vn1)
    λ = 0.9995 
    for i = 1:blksize:length(vn1)-blksize
        vn_best1 = compDoppler(vn1[i:i+blksize], yb[i:i+blksize])
        vn_best2 = compDoppler(vn2[i:i+blksize], yb[i:i+blksize])
        vn_best3 = compDoppler(vn3[i:i+blksize], yb[i:i+blksize])        
        if(length(vn_best2) > blksize)
            vn_best2 = vn_best2[1:blksize]
        end
        if(length(vn_best2) < blksize)
            vn_best2 = [vn_best2; zeros(ComplexF64,blksize-length(vn_best2))]
        end

        if(length(vn_best1) > blksize)
            vn_best1 = vn_best1[1:blksize]
        end
        if(length(vn_best1) < blksize)
            vn_best1 = [vn_best1; zeros(ComplexF64,blksize-length(vn_best1))]
        end

        if(length(vn_best3) > blksize)
            vn_best3 = vn_best3[1:blksize]
        end
        if(length(vn_best3) < blksize)
            vn_best3 = [vn_best3; zeros(ComplexF64,blksize-length(vn_best3))]
        end

        quantizer_op1, error, ff_filt1, SD, ffops, ff = chRlsTrain(vn_best1[1:blksize], yb[i:i+blksize], Q*12, λ, 1)
        Ber1 = computeBER(quantizer_op1, yb[i:i+blksize])
        quantizer_op2, error, ff_filt2, SD, ffops, ff = chRlsTrain(vn_best2[1:blksize], yb[i:i+blksize], Q*12, λ, 1)
        Ber2 = computeBER(quantizer_op2, yb[i:i+blksize])
        quantizer_op3, error, ff_filt3, SD, ffops, ff = chRlsTrain(vn_best3[1:blksize], yb[i:i+blksize], Q*12, λ, 1)
        Ber3 = computeBER(quantizer_op3, yb[i:i+blksize])
        println((i,Ber1, Ber2, Ber3))

        quantizer_op1, SD1 = DA(vn_best1[1:blksize], ff_filt1, 1)
        quantizer_op2, SD2 = DA(vn_best2[1:blksize], ff_filt2, 1)
        quantizer_op3, SD3 = DA(vn_best3[1:blksize], ff_filt3, 1)
        SDT = (0.3*SD1 + 0.5*SD2 + 0.2*SD3)
        atx=angle.(samples(yb[i:i+blksize]))
        tx_m=mean(atx)
        aSD = abs.(angle.(SDT))
        tx_zm = atx .- tx_m
        aSD_zm = aSD .- tx_m
        #println(length(aSD_zm))
        downbb, downtxbb = downsymboling(aSD_zm, tx_zm)
        ber = computeBER(downbb, downtxbb)
        println((i,ber))
        push!(e2,ber)
    end
    e2
end












y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
filename = "loc5-1.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
bbc, bb2 = prepsig(x[:,2])
events = Int[]
events = find_mseq(bbc, (yb), 0.7, 2000) #0.7, 3000 for 4-2 #0.5, 3000 for 5-1 # 0.7, 2000 for 5-2 # 0.7,2000
e_all = Float64[]
#events = events[2]
aaa = Float64[]
bbb = Float64[]
ccc = Int64[]
ddd = Float64[]
for i = 1:length(events)
    for j = 4:1:4
        e2 = Float64[]
        e3 = Float64[]
        temp = process(bbc, events[i], 1, j)
        temp2 = process3(events[i], 1, j)
        #println(length(temp))
        if(isassigned(temp))
            for v = 1:length(temp)
                push!(e2, temp[v])
                push!(e3, temp2[v])
            end
                if(!isempty(e2))
                    e2m = mean(e2)
                    e3m = mean(e3)
                    push!(aaa,e2m)
                    push!(bbb,j)
                    push!(ccc,i)
                    push!(ddd,e3m)
                end
        end
    end    
        # for k = 1:length(temp)
        #     push!(e_all, temp[k])
        # end
end
#plot(e2)
#plot!(repeat([mean(e2)],length(e2)))





