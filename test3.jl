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
ff_filt_len = 12*3
fb_filt_len = 12
SNR_dB = 30
snr = 10^(0.1*SNR_dB)
noise_var_1D = 0.5*2*1/(2*snr)
fade_chan = [repeat([0.9+0.9im],12); repeat([0.1-0.1im],12); repeat([-0.3+0.3im],12)]; 
fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan)
noise = sqrt(noise_var_1D)*randn(1, length(yb)+chan_len-1) + im*sqrt(noise_var_1D)*randn(1, length(yb)+chan_len-1)
chan_op = conv(fade_chan, samples(yb)) .+ noise'
chan_op = vec(chan_op)
chan_opre = resample(chan_op, 999/1000)
#chan_op = resample(chan_op, 999/1000)
quantizer_op, error, ff_filt ,fb_filt, SD1, ff1, fb1, aaa1= adaptRLSTrain2(chan_op[1:D], yb[1:D], ff_filt_len, fb_filt_len, 0.9, 1)
#quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op[1:D ], yb[1:D], ff_filt_len, fb_filt_len, 0.001, 1)
bits = Int64[]
[push!(bits, (sign(real(yb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(quantizer_op))]
ber1= 1 - (sum(bits)/length(bits))
println(ber1)
quantizer_op, error,ff_filt, SD, ffops = chRlsTrain(chan_op[1:D], yb[1:D], ff_filt_len, 0.99,1)
quantizer_op, SD2, fbops, ffops = adaptDD(chan_op[1:D], ff_filt , [0.04+0*im; 0.6+0*im] , 1)
bits = Int64[]
[push!(bits, (sign(real(yb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(quantizer_op))]
ber1= 1 - (sum(bits)/length(bits))
println(ber1)


# h= estimateD(chan_op, yb)/length(yb)
# h=h./norm(h,2)

# plot!(abs.(samples(h[1:8])))



plot(abs.(fade_chan))
ff_filt = ff_filt/norm(ff_filt,2)
plot!(abs.(reverse(ff_filt)))
 