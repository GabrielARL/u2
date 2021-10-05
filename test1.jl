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

#include("gabUtil.jl")
wavsize(filename) = wavread(filename; format="size")
fs = 32000
bpf = fir(256, 3800.0, 6200.0; fs=fs)


# x = signal(filename; start=1, nsamples=min(blksize, stop-1+1))

# bits = Int64[]
# for i = 12:12:length(t)
#     push!(bits,convert.(Int64,mean(t[i-11:i])))
# end
# rbb_1 = ComplexF64[]
# for i = 12:12:length(rbb)
#     push!(rbb_1,convert.(ComplexF64,mean(rbb[i-11:i])))
# end

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
 #this is the adaptive LMS DFE  
 function adaptLMSTrain(s, bb, N_ff, N_fb, δ)
  ff_filt_len = N_ff
  fb_filt_len = N_fb
  ff_filt = Array{ComplexF64}(undef,1,ff_filt_len) 
  fb_filt = Array{ComplexF64}(undef,1,fb_filt_len)  
  ff_filt_ip = Array{ComplexF64}(undef,1,ff_filt_len)
  fb_filt_ip = Array{ComplexF64}(undef,1,fb_filt_len)

  ff_filt = zeros(ComplexF64,1,ff_filt_len)
  fb_filt = zeros(ComplexF64,1,fb_filt_len)
  ff_filt_ip = zeros(ComplexF64,1,ff_filt_len)
  fb_filt_ip = zeros(ComplexF64,1,fb_filt_len)

  fb_filt_op = 0
  quantizer_op = Array{ComplexF64}(undef,1,(length(bb)-ff_filt_len+1))
  err = Array{ComplexF64}(undef,1,(length(bb)-ff_filt_len+1))
  for i1=1:length(bb)-ff_filt_len+1 
    ff_filt_ip[2:end]=ff_filt_ip[1:end-1];
    ff_filt_ip[1] = s[i1];
    ff_filt_op = vec(ff_filt)'*vec(ff_filt_ip)
    ff_and_fb = ff_filt_op-fb_filt_op; 
    err[i1] = ff_and_fb - bb[i1]  
    #println((i1, fb_filt_op))
    #temp = imag(ff_and_fb) < 0.0 ? 0 : 1
    quantizer_op[i1] = sign(real(ff_and_fb))
    ff_filt=ff_filt-δ*err[i1]*ff_filt_ip;
    fb_filt=fb_filt+δ*err[i1]*fb_filt_ip;
    fb_filt_ip[2:end]=fb_filt[1:end-1];
    fb_filt_ip[1] = quantizer_op[i1];
    fb_filt_op = vec(fb_filt)'*vec(fb_filt_ip);
    if(any(isinf,fb_filt_op))
      #println((i1, err[i1], fb_filt_op, any(isnan,fb_filt), any(isnan,ff_filt_ip)))
      break
    end 
  end
  collect(vec(quantizer_op)), ff_filt, fb_filt, collect(vec(err))
end

function adaptRLSTrain(s, bb, N_ff, N_fb, λ, Δ)
  ff_filt_len = N_ff
  fb_filt_len = N_fb
  ff_filt = Array{ComplexF64}(undef,ff_filt_len) 
  fb_filt = Array{ComplexF64}(undef,fb_filt_len)  
  ff_filt_ip = Array{ComplexF64}(undef,ff_filt_len)
  fb_filt_ip = Array{ComplexF64}(undef,fb_filt_len)

  ff_filt = zeros(ComplexF64,ff_filt_len)
  fb_filt = zeros(ComplexF64,fb_filt_len)
  ff_filt_ip = zeros(ComplexF64,ff_filt_len)
  fb_filt_ip = zeros(ComplexF64,fb_filt_len)

  fb_filt_op = 0
  quantizer_op = Array{ComplexF64}(undef,(length(s)-ff_filt_len-fb_filt_len))
  SD = Array{ComplexF64}(undef,(length(s)-ff_filt_len-fb_filt_len))
  P1 = Δ * I(N_ff)
  P2 = Δ * I(N_fb)

  err = Array{ComplexF64}(undef,(length(s)-ff_filt_len-fb_filt_len))
  #println(length(s)-ff_filt_len-fb_filt_len)
  for i1=1:length(s)-ff_filt_len-fb_filt_len
    ff_filt_ip[2:end]=ff_filt_ip[1:end-1]
    ff_filt_ip[1] = s[i1]
    ff_filt_op = vec(ff_filt)'*vec(ff_filt_ip)
    ff_and_fb = ff_filt_op-fb_filt_op; 
    err[i1] =  -(bb[i1] - ff_and_fb) 
    #println((i1, fb_filt_op))
    #temp = imag(ff_and_fb) < 0.0 ? 0 : 1
    quantizer_op[i1] = sign(real(ff_and_fb))
    SD[i1] = ff_and_fb
    #println((i1, ff_and_fb))
    u1 = ff_filt_ip
    nume1 = P1*u1
    denom1 = λ + u1'*nume1
    K1 = nume1/denom1
    PPrime1 = K1*u1'*P1
    PPrime1 = collect(PPrime1)
    P1=(P1.-PPrime1)/λ
    
    u2 = vec(fb_filt_ip)
    nume2 = P2*u2
    denom2 = λ + u2'*nume2
    K2 = nume2/denom2
    PPrime2 = K2*u2'*P2
    PPrime2 = collect(PPrime2)
    if( N_fb == 1 )
      P2= (only(P2)-only(PPrime2)) /λ
    else
      P2=(P2-(PPrime2))/λ
    end
    ff_filt=ff_filt-K1*err[i1];
    if( N_fb == 1)
      fb_filt=only(fb_filt)+only(K2*err[i1]);
      fb_filt_ip[1] = quantizer_op[i1];
      fb_filt_op = only(fb_filt * fb_filt_ip)
    else
      fb_filt=fb_filt+K2*err[i1];
      fb_filt_ip[2:end]=fb_filt[1:end-1];
      fb_filt_ip[1] = quantizer_op[i1];
      fb_filt_op = vec(fb_filt)'*vec(fb_filt_ip);
    end
    #fb_filt_op = vec(fb_filt)'*vec(fb_filt_ip);
  end
  quantizer_op[1:end], ff_filt, fb_filt, err, SD
end


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
#display(plot(abs.(mfilter(yb[1:3000],bb))))
events = Int[]
events=find_mseq(bb, real.(yb))
T = length(yb)

containers = Array{ComplexF64}[]
containers1 = Array{ComplexF64}[]
empty!(containers)
empty!(containers1)
for i ∈ eachindex(events)
    rbb = bb[max(events[i], 1) : events[i] + T  ] 
    push!(containers, rbb)
    #eqlz, w, wb, e, SD=adaptRLSTrain(rbb[1:3000], yb[1:3000], 10, 5, 0.999, 0.01)
end



# bits = Int64[]
# for i = 12:12:length(yb)
#     push!(bits,convert.(Int64,mean(yb[i-11:i])))
# end

tx_bb = repeat(gmseq(12),1)
tx_bb = round.(tx_bb - real(tx_bb))
tx_bb = -im*tx_bb

for p ∈ eachindex(containers)
rbb_1 = ComplexF64[]

for i = 12:12:length(rbb)
    push!(rbb_1,convert.(ComplexF64,mean(containers[p][max(1,i-11):i])))
end

push!(containers1,rbb_1)

end


# channel 2 is the best, loc 4-2 can see, loc5-1 loc5-2, loc6

# h=MMSE_eqlz(rbb_1, tx_bb, 10)
# sd=sfilt(h, 1, rbb_1)
# sd1 = -sd[kk:end]
eqlz, w, wb, e, SD=adaptRLSTrain(rbb_1, tx_bb, 3, 3, 0.99, 0.001)
eqlz, w, wb, e = adaptLMSTrain(rbb_1, tx_bb, 10, 1, 0.001)
sd1 = SD[1:end]
bits = Int[]
#(plot(real.(tx_bb[10^3:2*10^3])))
display(plot!(sign.(real.(sd1[10^3:2*10^3]))))
#display(plot!(abs.(e)))
[push!(bits, sign(real(tx_bb[i]))==sign(real(sd1[i])) ? 1 : 0)  for i ∈ 1:min(length(sd1),length(tx_bb))] 
ber1= 1 - (sum(bits)/length(bits))
if(ber1 < 0.47)
println(ber1)
end





#sd1 = -sd[1:end-13]



