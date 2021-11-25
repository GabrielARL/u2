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
ff_filt_len = 111
fb_filt_len = 20
SNR_dB = 30
snr = 10^(0.1*SNR_dB)
noise_var_1D = 0.5*2*1/(2*snr)
fade_chan = [0.1+0.1im 0.1+0.1im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im 0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im 1.2-3.2im 1.1-1.1im 0.2-0.2im 3.0 - 1.2im 1.4+1.2im 1.1+1.3im]'; 
#fade_chan = fade_chan*0.01

fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);
noise = sqrt(noise_var_1D)*randn(1, length(yb)+chan_len-1) + im*sqrt(noise_var_1D)*randn(1, length(yb)+chan_len-1)
chan_op = conv(fade_chan, samples(yb)) .+ noise'
chan_op = vec(chan_op)
#chan_opre = resample(chan_op, 999/1000)
quantizer_op, error, ff_filt ,fb_filt, SD1, ff1, fb1, aaa1= adaptRLSTrain2(chan_op, samples(yb), ff_filt_len, fb_filt_len, 0.9999, 1)
#quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op[1:D ], yb[1:D], ff_filt_len, fb_filt_len, 0.001, 1)
bits = Int64[]
[push!(bits, (sign(real(yb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(quantizer_op))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))
quantizer_op, SD2, fbops, ffops = adaptDD(chan_op, ff_filt ,fb_filt, 1)
bits = Int64[]
[push!(bits, (sign(real(yb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(quantizer_op))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))

# txbb=samples(mseq(12))
# chdown=downsymbol(chan_opre[1:end])
# quantizer_op, error, ff_filt ,fb_filt, SD1, ff1, fb1, aaa1= adaptRLSTrain2(chdown, txbb, ff_filt_len, fb_filt_len, 0.9999, 1)
# #quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op[1:D ], yb[1:D], ff_filt_len, fb_filt_len, 0.001, 1)
# bits = Int64[]
# [push!(bits, (sign(real(txbb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(txbb),length(quantizer_op))]
# ber1= 1 - (sum(bits)/length(bits))
# println((ber1))
# quantizer_op, SD2, fbops, ffops = adaptDD(chdown, ff_filt ,fb_filt, 1)
# bits = Int64[]
# [push!(bits, (sign(real(txbb[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(txbb),length(quantizer_op))]
# ber1= 1 - (sum(bits)/length(bits))
# println((ber1))

vbb=resample(chan_op, 1//12)
xbb=mseq(12)
h= estimateD(vbb, xbb)/length(xbb)
h=h./norm(h,2)
f = resample(vec(fade_chan), 1//12)
f = f/length(f)
f = f/norm(f,2)
