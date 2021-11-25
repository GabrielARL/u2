using SignalAnalysis
using Statistics
using LinearAlgebra
using Random
using Plots
using DSP

include("rls.jl")

xbb = repeat(mseq(12),10)

ff_filt_len = 3
fb_filt_len = 2
SNR_dB = 30
snr = 10^(0.1*SNR_dB)
noise_var_1D = 0.5*2*1/(2*snr)
fade_chan = [0.1+0.1im 0.1+0.1im 0.9+0.9im 0+0im 0+0im]'

fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);
noise = sqrt(noise_var_1D)*randn(1, length(xbb)+chan_len-1) + im*sqrt(noise_var_1D)*randn(1, length(xbb)+chan_len-1)
chan_op = conv(fade_chan, samples(xbb)) .+ noise'
chan_op = vec(chan_op)

h= estimateD(chan_op, xbb)/length(xbb)
h=h./norm(h,2)
plot(abs.(fade_chan))
plot!(abs.(samples(h[1:6])))
h = MMSE_eqlz( chan_op, samples(xbb), 6)
h = h/norm(h, 2)
plot!(abs.(samples(h[1:6])))

quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op, xbb, 6, 6, 0.001, 1)
bits = Int64[]
[push!(bits, (sign(real(xbb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(xbb),length(quantizer_op))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))
ff_filt = ff_filt/norm(ff_filt,2)
plot!(abs.(reverse(ff_filt)))
fb_filt = fb_filt/norm(fb_filt,2)
plot!(abs.(reverse(fb_filt)))