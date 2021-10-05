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

#include("utils.jl")


filename = "loc6.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y.*cw(1000, length(mseq(12))*12/6000, 6000.0)
#y = repeat(mseq(12); inner=12)
d1 = sfiltfilt(bpf, x[:,1])
d2 = sfiltfilt(bpf, x[:,2])
d3 = sfiltfilt(bpf, x[:,3])
display(specgram(d2;fs=32000,nfft=4096, noverlap=2048))
bb = downconvert(sresample(d2, 9//8), 6, 6000)
#plot(abs.(mfilter(y[1:3000],bb)))
events = Int[]
events=find_mseq(bb, y)

T = length(y)
for i ∈ events
    rbb = bb[events[i]:events[i] + T + 100] 
    #h = MMSE_eqlz(rbb, y, 10)
    #println(i)
    
end








# cd(dirname(Pkg.project().path))
# filename = "loc4-1.wav" 
# nsamples, _ = wavsize(filename)

# stop = nsamples
# start = 1
 wavsize(filename) = wavread(filename; format="size")

# z = Float64[]
# #fs = 0.0
fs = 32000
bpf = fir(256, 3800.0, 6200.0; fs=fs)


# x = signal(filename; start=1, nsamples=min(blksize, stop-1+1))
# y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
# y = repeat(mseq(12); inner=12)
# bits = Int64[]
# for i = 12:12:length(t)
#     push!(bits,convert.(Int64,mean(t[i-11:i])))
# end
# rbb_1 = ComplexF64[]
# for i = 12:12:length(rbb)
#     push!(rbb_1,convert.(ComplexF64,mean(rbb[i-11:i])))
# end
function find_mseq(bb, tx)
    λ = 0.9999          # exponential averaging factor for threshold
    β = 5.0             # threshold is β × average
    xmin = 0.10       # minimum threshold
    gap = 3500         # minimum gap between detections
    pwidth = 25        # peak width (to look for maxima)
    x = abs.(mfilter(tx[1:3000], bb))
    μ = x[1]
    j = 0
    events = Int[]
    fs = 6000
    for i = 1:length(x)
      μ = λ * μ + (1 - λ) * x[i]
      if x[i] ≥ max(xmin, β * μ) && i ≥ j
        push!(events, argmax(x[i:min(length(x),i+pwidth)]) + i - 1)
        j = i + gap
      end
    end
    invalid = ones(Bool, size(events))
    for i = 2:length(events)
      dt = (events[i] - events[i-1]) / fs
      n = dt/19                                       # 19s gap between LFMs in a group
      if n > 1.0 && n < 2.0 && isapprox(n, round(n); atol=0.2)  # valid -> gap is right
        invalid[i] = false
        invalid[i-1] = false
      end
    end
    println(invalid)
    deleteat!(events, findall(invalid))
    events
  end 

  function MMSE_eqlz(s, MSEQ, M)
    ϕff = xcorr(samples(s), samples(s))
    ϕfd = xcorr(MSEQ, samples(s))
    N = length(MSEQ)
    𝑅ff = ϕff[N:N+M-1]
    𝑅 = Toeplitz(𝑅ff, vec(𝑅ff))
    P = ϕfd[N:N+M-1]
    h = inv(collect(𝑅)) * (P)
    h
  end









# channel 2 is the best, loc 4-2 can see, loc5-1 loc5-2, loc6

