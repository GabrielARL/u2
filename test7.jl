using SignalAnalysis
using Statistics
using LinearAlgebra
using Random
using Plots
using DSP

include("rls.jl")

y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
D = 21000
ff_filt_len = 12
fb_filt_len = 12
SNR_dB = 30
snr = 10^(0.1*SNR_dB)
noise_var_1D = 0.5*2*1/(2*snr)
fade_chan = [0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im  0.9+0.9im 0.1 +0.2im]'; 
#fade_chan = real(fade_chan)
fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);
noise = sqrt(noise_var_1D)*randn(1, length(yb)+chan_len-1) + im*sqrt(noise_var_1D)*randn(1, length(yb)+chan_len-1)
chan_op = conv(fade_chan, samples(yb)) .+ noise'
chan_op = vec(chan_op)
chan_opre = resample(chan_op, 999/1000)
h = MMSE_eqlz( chan_op, samples(yb), 6)
h = h/norm(h, 2)
sd = sfilt(h, 1, chan_op)
bits = Int64[]
[push!(bits, (sign(real(yb[i]))==sign(real(sd[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(sd))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))
quantizer_op, error,ff_filt, SD, ffops = chRlsTrain(chan_op, yb, 6, 0.99,1)
bits = Int64[]
ff_filt = ff_filt/norm(ff_filt,2)
[push!(bits, (sign(real(yb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(quantizer_op))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))
plot(abs.(fade_chan))
plot!(abs.(h))
plot!(abs.(ff_filt))

