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
wavsize(filename) = wavread(filename; format="size")
fs = 32000
bpf = fir(256, 3800.0, 6200.0; fs=fs)

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

    println(events)

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

y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)

filename = "loc4-2.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
#idx = findall(x-> abs(x) > 0.02, x[:,2])
#x[idx,2] .= 0 
d1 = sfiltfilt(bpf, x[:,1])
d2 = sfiltfilt(bpf, x[:,2])
d3 = sfiltfilt(bpf, x[:,3])
#display(specgram(d2;fs=32000,nfft=4096, noverlap=2048))
bb=downconvert(sresample(d2, 9//8), 6, 6000)
bb = bb .* cw(1000.0, length(bb)/6000, 6000.0)

bb2=downconvert(sresample(d1, 9//8), 6, 6000)
bb2 = bb2 .* cw(1000.0, length(bb2)/6000, 6000.0)

bb3=downconvert(sresample(d3, 9//8), 6, 6000)
bb3 = bb3 .* cw(1000.0, length(bb3)/6000, 6000.0)

bbc=downconvert(sresample(x[:,2], 9//8), 6, 6000)
bbc = bbc .* cw(1000.0, length(bbc)/6000, 6000.0)

events = Int[]
events = find_mseq(bb, real.(yb))
T = length(yb)

txbb=convert.(Int64, real.(mseq(12)));


rbb = bbc[events[1]:events[1]+T]
rbb = resample(rbb,1000/999)
rbb = rbb./norm(rbb,2)
c = MMSE_eqlz(rbb, samples(yb), 12*3)
h = (1/c) 
h = h/norm(h,2)
h = h'

vn = zeros(ComplexF64, length(rbb), 1)
for i =length(h):length(rbb)
vn[i] = rbb[i] - vec(h[1:end])'*reverse(rbb[i-length(h)+1:1:i])
end

vn1 = resample(vn[1:min(length(yb), length(vn))],(1000)/1000)
#quantizer_op, error, ff_filt, SD, ffops, ff = chRlsTrain(vn1[1:3000], yb[1:length(vn1)], 3*12, 0.99, 1)
# bits = Int64[]
# [push!(bits, (sign(real(yb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(quantizer_op))]
# ber1= 1 - (sum(bits)/length(bits))
#println((ber1))
# ff_filt = ff_filt/norm(ff_filt,2)
# println(ff_filt)
#quantizer_op = -quantizer_op

e = Float64[]
e2 = Float64[]
SDs = Float64[]
atxDs = Float64[]

downbb = Int64[]
downtxbb = Int64[]

for i = 1:1200:length(vn)-1200
  a = Float64[]
  for j = 1:1:200
    vn2 = resample(vn1[i:i+1200], (j+900)/1000)
    cr=mfilter(yb[i:i+1200], vn2)
    max = maximum(abs.(cr))
    push!(a, max)
  end
fac=findmax(a)[2]
#println(fac)
vn_best = resample(vn1[i:i+1200], (fac+900)/1000)
if(length(vn_best) > 1200)
  vn_best = vn_best[1:1200]
end
if(length(vn_best) < 1200)
  vn_best = [vn_best; zeros(ComplexF64,1200-length(vn_best))]
end
quantizer_op, error, ff_filt, SD, ffops, ff = chRlsTrain(vn_best[1:1200], yb[i:i+1200], 4*12, 0.9995, 1)
quantizer_op, SD2 = DA(vn_best[1:1200], ff_filt, 1)
atx=angle.(samples(yb[i:i+1200]))
tx_m=mean(atx)
aSD = abs.(angle.(SD2))
tx_zm = atx .- tx_m
aSD_zm = aSD .- tx_m
bits = Int64[]
[push!(bits, (sign((tx_zm[i]))==sign((aSD_zm[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(SD2))]
ber1= 1 - (sum(bits)/length(bits))
#println((i,ber1))
push!(e,ber1)



for h = 1:length(SD2)
  push!(SDs, aSD_zm[h])
  push!(atxDs, tx_zm[h])
  if(h%12==0)
    temp1 = sign.(aSD_zm[h:-1:(h-12+1)])
    n=count(>=(0.0), temp1)
    if (n >= 6)
      bit = 1
    else 
      bit = -1
    end
    push!(downbb,bit)
    temp1 = sign.(tx_zm[h:-1:(h-12+1)])
    n=count(>=(0.0), temp1)
    println(n)
    if (n >= 6)
      bit = 1
    else 
      bit = -1
    end
    push!(downtxbb,bit)
  end
end

bits = Int64[]
[push!(bits, (sign((downbb[i]))==sign((downtxbb[i])))  ? 1 : 0)  for i ∈ 1:min(length(downbb),length(downtxbb))]
ber2= 1 - (sum(bits)/length(bits))
println((i,ber1,ber2))
push!(e2,ber2)
end













# push!(a,i)
# push!(b,j)
# push!(c,k)

# txbb=convert.(Int64,real.(samples(round.(yb))))
# delay = 1
# rbb = bb[max(events[1],1)+delay:events[1]+T+delay] 
# rbb_d = resample(rbb, 1//12)
# tx=mseq(12)




