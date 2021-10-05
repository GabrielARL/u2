using ToeplitzMatrices: size
using StatsPlots: length, repeat, mean
include("gabUtil.jl")
#plotlyjs()
datasets  = 4 #[1:2]
data, timings = getTimings(datasets)
dataDf, timings = packSignals(data,timings)
bb = repeat(gmseq(8), 1)
bb = round.(bb - real(bb))
bb = -im*bb
bb = repeat(bb,5)
ch = 1
chs = 12
tr_len = convert(Int64,round(1*length(bb)))
temp = ComplexF64[]
container = [[1.0+1.0im, 1.0+1.0im]]
containerSD = [[1.0+1.0im, 1.0+1.0im]]
bers = Float64[]
empty!(containerSD)

for k = 1:chs
for i = 1:maximum(timings[:group])
    for j = 1:k
        if(!isnothing(dataDf[i,j,1]))
            sig = dataDf[i,j,:]
            sig1=convert.(Float64, sig)
            x = sresample(sig1, 44100/fs)
            x = x[1:end-1]
            x = samples(x)
            #println((i,j))
            push!(container, x)
            if (i == 1 && j == 1)
                deleteat!(container, 1)
            end
        end
    end
        if (size(container)[1] >= 1) 
            if(size(container[1])[1] == 2)
                deleteat!(container, 1)
            end
            for i = 1:size(container)[1]
                eqlz, w, wb, e, SD=adaptRLSTrain(container[i][1:tr_len], bb[1:tr_len], 4, 4, 0.85, 0.05)
                #eqlz, SD = adaptRLSDD(container[i][1:tr_len], w ,wb)
                push!(containerSD, SD[2:end])
            end
        if(size(containerSD)[1] > 0)
        eqlz = sum(containerSD[1:size(containerSD)[1]])/size(containerSD)[1]
        empty!(containerSD)
        dec_a = sign.(real.(eqlz)) 
        dec = convert.(Int64, dec_a)
        dec = dec[1:end]
        bits = Int[]
        data1 = real.(bb[tr_len:end])
        [push!(bits, bb[i]==dec[i] ? 1 : 0)  for i ∈ 1:length(dec)] 
        ber=1-(sum(bits)/length(bits))
        #println(ber)
        push!(bers,ber)
        end
        end
    #println(size(container))
    empty!(container)
    
end
println((mean(bers),k))
empty!(bers)
GC.gc()
end
 