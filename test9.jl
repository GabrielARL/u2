using SignalAnalysis
using Statistics
using LinearAlgebra
using Random
using Plots
using DSP

include("rls.jl")
xb=convert.(ComplexF64, repeat(mseq(12),1))
#xb = convert.(Float64, real.(xb))
SNR_dB = 30
snr = 10^(0.1*SNR_dB)
noise_var_1D = 0.5*2*1/(2*snr)
fade_chan = [0.8+0.2*im 0.1+0.1*im 0.1+0.3*im]';
fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);
noise = sqrt(noise_var_1D)*randn(1, length(xb)+chan_len-1)
chan_op = conv(fade_chan, xb) 
chan_op = vec(chan_op)
#chan_op = resample(chan_op, 999/1000)
h= estimateD(chan_op, xb)/length(xb)
plot(fade_chan)
plot!(real(samples(h[1:3])))
vn = zeros(ComplexF64, length(chan_op),1)
vn1 = zeros(ComplexF64, length(chan_op),1)
for i = 3:length(chan_op)
vn[i] = chan_op[i] - (h[1]*chan_op[i] + h[2]*chan_op[i-1] + h[3]*chan_op[i-2])
vn1[i] = chan_op[i] - (samples(h[1:3]))'*reverse(chan_op[i-2:1:i])
end
quantizer_op, error, ff_filt, SD, ffops = chRlsTrain(vn, xb, 3, 0.99, 1)
xb = xb[1:end]
bits = Int64[]
[push!(bits, (sign(real(xb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(xb),length(quantizer_op))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))