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

filename = "loc6.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
#idx = findall(x-> abs(x) > 0.02, x[:,2])
#x[idx,2] .= 0 
d1 = sfiltfilt(bpf, x[:,1])
d2 = sfiltfilt(bpf, x[:,2])
d3 = sfiltfilt(bpf, x[:,3])
display(specgram(d1;fs=32000,nfft=4096, noverlap=2048))
display(specgram(d2;fs=32000,nfft=4096, noverlap=2048))
display(specgram(d3;fs=32000,nfft=4096, noverlap=2048))
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

containers = Array{ComplexF64}[]
containers1 = Array{ComplexF64}[]
empty!(containers)
empty!(containers1)
txbb=convert.(Int64,real.(samples(round.(yb))))
Q = 7;
ff_filt_len = Q*12
fb_filt_len = Q*12+1
bb = bbc
delay = 6
rbb = bb[max(events[1], 1)+delay : events[1]+T+delay] 
rbb2 = bb2[max(events[1], 1)+delay : events[1]+T+delay]
rbb3 = bb3[max(events[1], 1)+delay : events[1]+T+delay]
gain = 1
rbb_norm = rbb*norm(txbb,2)/norm(rbb,2)*gain
rbb_norm = sresample(rbb_norm, 1000/990)*gain
rbb_norm2 = rbb2*norm(txbb,2)/norm(rbb2,2)*gain
rbb_norm2 = sresample(rbb_norm2, 1000/990)*gain
rbb_norm3 = rbb3*norm(txbb,2)/norm(rbb3,2)*gain
rbb_norm3 = sresample(rbb_norm3, 1000/990)*gain
rbb_normD = rbb_norm[1:6000]
rbb_normT = rbb_norm[1:6000]
rbb_normD2 = rbb_norm2[1:6000]
rbb_normT2 = rbb_norm2[1:6000]
rbb_normD3 = rbb_norm3[1:6000]
rbb_normT3 = rbb_norm3[1:6000]

quantizer_op, error,ff_filt, fb_filt, SD = adaptRLSTrain(samples(rbb_normT), txbb[1:length(rbb_normT)], ff_filt_len, fb_filt_len, 0.999, 1)
quantizer_op1, SD1, fbops, ffops=adaptDD(samples(rbb_normD), ff_filt*0.88 ,fb_filt*0.1, 1)
bits = Int[]
bits = Int[];[push!(bits, sign(real(txbb[i]))==sign(real(quantizer_op1[i])) ? 1 : 0)  for i ∈ 1:min(length(quantizer_op1),length(txbb))] 
ber1= 1 - (sum(bits)/length(bits))
println(ber1)

quantizer_op, error,ff_filt, fb_filt, SD = adaptRLSTrain(samples(rbb_normT2), txbb[1:length(rbb_normT)], ff_filt_len, fb_filt_len, 0.999, 1)

quantizer_op2, SD2, fbops, ffops=adaptDD(samples(rbb_normD2), ff_filt*0.99 ,fb_filt*0.05, 1)
bits = Int[]
bits = Int[];[push!(bits, sign(real(txbb[i]))==sign(real(quantizer_op2[i])) ? 1 : 0)  for i ∈ 1:min(length(quantizer_op2),length(txbb))] 
ber2= 1 - (sum(bits)/length(bits))
println((ber2))

quantizer_op, error,ff_filt, fb_filt, SD = adaptRLSTrain(samples(rbb_normT3), txbb[1:length(rbb_normT)], ff_filt_len, fb_filt_len, 0.999, 1)
quantizer_op3, SD3, fbops, ffops=adaptDD(samples(rbb_normD3), ff_filt*0.98 ,fb_filt*0.07, 1)
bits = Int[];[push!(bits, sign(real(txbb[i]))==sign(real(quantizer_op3[i])) ? 1 : 0)  for i ∈ 1:min(length(quantizer_op3),length(txbb))] 
ber3= 1 - (sum(bits)/length(bits))
println((ber3))

spatial_comb_decision = Int[];[push!(spatial_comb_decision, sign(round(quantizer_op3[i]+quantizer_op2[i]+quantizer_op1[i]) ))  for i ∈ 1:min(length(quantizer_op1),length(txbb))] 
bits = Int[];[push!(bits, sign(real(txbb[i]))==sign(real(spatial_comb_decision[i])) ? 1 : 0)  for i ∈ 1:min(length(quantizer_op2),length(txbb))]
ber4= 1 - (sum(bits)/length(bits))
println((ber4))
ttt=resample(spatial_comb_decision, 1//12)
SDT= SD1*norm(SD1,2)+ SD2*norm(SD2,2) + SD3*norm(SD3,2)
SDTT=sign.(real.(SDT))
bits = Int[];[push!(bits, sign(real(txbb[i]))==sign(real(SDTT[i])) ? 1 : 0)  for i ∈ 1:min(length(quantizer_op2),length(txbb))]
ber5= 1 - (sum(bits)/length(bits))
println((ber5))

# for j = 0:12
# println((0,0))
# b=sign.(real.(sresample(SD[4:end],1//12)))
# b = b[1+j:end]
# a=mseq(12)[1:length(b)]
# bits = Int[];[push!(bits, sign(real(a[i]))==sign(real(b[i])) ? 1 : 0)  for i ∈ 1:min(length(b),length(a))] 
# ber1= 1 - (sum(bits)/length(bits))
# println((j, ber1))
# end


# bits = Int64[]
# for i = 12:12:length(yb)
#     push!(bits,convert.(Int64,mean(yb[i-11:i])))
# end

# tx_bb = repeat(gmseq(12),1)
# tx_bb = round.(tx_bb - real(tx_bb))
# tx_bb = -im*tx_bb

# for p ∈ eachindex(containers)
# rbb_1 = ComplexF64[]

# for i = 12:12:length(rbb)
#     push!(rbb_1,convert.(ComplexF64,mean(containers[p][max(1,i-11):i])))
# end

# push!(containers1,rbb_1)

# end


# channel 2 is the best, loc 4-2 can see, loc5-1 loc5-2, loc6

# h=MMSE_eqlz(rbb_1, tx_bb, 10)
# sd=sfilt(h, 1, rbb_1)
# sd1 = -sd[kk:end]
# eqlz, w, wb, e, SD=adaptRLSTrain(rbb_1, tx_bb, 3, 3, 0.99, 0.001)
# eqlz, w, wb, e = adaptLMSTrain(rbb_1, tx_bb, 10, 1, 0.001)
# sd1 = SD[1:end]
# bits = Int[]
# #(plot(real.(tx_bb[10^3:2*10^3])))
# display(plot!(sign.(real.(sd1[10^3:2*10^3]))))
# #display(plot!(abs.(e)))
# [push!(bits, sign(real(tx_bb[i]))==sign(real(sd1[i])) ? 1 : 0)  for i ∈ 1:min(length(sd1),length(tx_bb))] 
# ber1= 1 - (sum(bits)/length(bits))
# if(ber1 < 0.47)
# println(ber1)
# end





#sd1 = -sd[1:end-1]