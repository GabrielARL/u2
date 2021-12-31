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





