using SignalAnalysis
using Statistics
using LinearAlgebra
using Random
using Plots
using DSP

include("rls.jl")

xb=convert.(ComplexF64, repeat(mseq(12),2))
SNR_dB = 30
snr = 10^(0.1*SNR_dB)
noise_var_1D = 0.5*2*1/(2*snr)
fade_chan = [0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im  0.9+0.9im 0.1 +0.2im]'; 
fade_chan = repeat(fade_chan,3)
#fade_chan = real(fade_chan)
fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);
noise = sqrt(noise_var_1D)*randn(1, length(xb)+chan_len-1) + im*sqrt(noise_var_1D)*randn(1, length(xb)+chan_len-1)
chan_op = conv(fade_chan, samples(xb)) .+ noise'
chan_op = vec(chan_op)
chan_opre = resample(chan_op, 999/1000)
h = MMSE_eqlz( chan_op, samples(xb), 21)
h = h/norm(h, 2)
sd = sfilt(h, 1, chan_op)
bits = Int64[]
[push!(bits, (sign(real(xb[i]))==sign(real(sd[i])))  ? 1 : 0)  for i ∈ 1:min(length(xb),length(sd))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))
quantizer_op, error,ff_filt, SD, ffops = chRlsTrain(chan_op, yb, 21, 0.99,1)
bits = Int64[]
ff_filt = ff_filt/norm(ff_filt,2)
plot(abs.(fade_chan))
plot!(abs.(h))
plot!(abs.(reverse(ff_filt)))



#quantizer_op, error, ff_filt ,fb_filt, SD1, ff1, fb1, aaa1= adaptRLSTrain2(chan_op, samples(xb), 21, 5, 0.99, 1)
#quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op[1:D ], yb[1:D], ff_filt_len, fb_filt_len, 0.001, 1)
# bits = Int64[]
# [push!(bits, (sign(real(xb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(xb),length(quantizer_op))]
# ber1= 1 - (sum(bits)/length(bits))
# println((ber1))
# ff_filt = ff_filt/norm(ff_filt,2)
#plot!(abs.(reverse(ff_filt)))


# quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op, xb, 21, 7, 0.001, 1)
# bits = Int64[]
# [push!(bits, (sign(real(xb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(xb),length(quantizer_op))]
# ber1= 1 - (sum(bits)/length(bits))
# println((ber1))
# ff_filt = ff_filt/norm(ff_filt,2)
#plot!(abs.(reverse(ff_filt)))