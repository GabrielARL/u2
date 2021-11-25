using Statistics
using LinearAlgebra
using Random
using SignalAnalysis
using DSP
using UnderwaterAcoustics
using Plots
Random.seed!(1234)

include("rls.jl")

training_len = 1*10^5
SNR_dB = 30
ff_filt_len = 30
fb_filt_len = 20
data_len = 1*10^5

snr = 10^(0.1*SNR_dB)
noise_var_1D = 0.5*2*1/(2*snr)
δ = 0.001
training_a = rand([0 1],1,2*training_len)
training_seq = 1 .- (2*training_a[1:2:end]) + 1im*(1 .- 2*training_a[2:2:end])

fade_chan = [0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im]'; 
fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);

noise = sqrt(noise_var_1D)*randn(1, training_len+chan_len-1) + im*sqrt(noise_var_1D)*randn(1, training_len+chan_len-1)
chan_op = conv(fade_chan, training_seq) .+ noise'

quantizer_op, ff_filt1, fb_filt1, err = adaptLMSTrain(chan_op, training_seq, ff_filt_len, fb_filt_len, 0.001, 2)
#quantizer_op, error, ff_filt1 ,fb_filt1,SD1 = adaptRLSTrain(chan_op, training_seq, ff_filt_len, fb_filt_len, 0.99, 2)
bits = Int64[]
[push!(bits, (sign(real(training_seq[i]))==sign(real(quantizer_op[i]))) & (sign(imag(training_seq[i]))==sign(imag(quantizer_op[i]))) ? 1 : 0)  for i ∈ 1:min(length(training_seq),length(quantizer_op))] 
ber1= 1 - (sum(bits)/length(bits))
println(ber1)

quantizer_op,SD2, fbops, ffops = adaptDD(chan_op, ff_filt1 ,fb_filt1, 2)
bits = Int64[]
[push!(bits, (sign(real(training_seq[i]))==sign(real(quantizer_op[i]))) & (sign(imag(training_seq[i]))==sign(imag(quantizer_op[i]))) ? 1 : 0)  for i ∈ 1:min(length(training_seq),length(quantizer_op))] 
ber1= 1 - (sum(bits)/length(bits))
println(ber1)

training_a = rand([0 1],1,training_len)
training_seq = 1 .- (2*training_a)' 
fade_chan = [0.9+0.9im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.1+0.1im 0.9+0.9im 1.2-0.3im]'; 
fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);

noise = sqrt(noise_var_1D)*randn(1, training_len+chan_len-1) + im*sqrt(noise_var_1D)*randn(1, training_len+chan_len-1)
chan_op = conv(fade_chan, training_seq) .+ noise'
quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op, training_seq, ff_filt_len, fb_filt_len, 0.005, 1)
#quantizer_op, error, ff_filt ,fb_filt,SD1 = adaptRLSTrain(chan_op, training_seq, ff_filt_len, fb_filt_len, 0.99, 1)
bits = Int64[]
[push!(bits, (sign(real(training_seq[i]))==sign(real(quantizer_op[i])))  ? 1 : 0)  for i ∈ 1:min(length(training_seq),length(quantizer_op))] 
ber1= 1 - (sum(bits)/length(bits))
println(ber1)
ff_filt=reverse(ff_filt)
quantizer_op, SD2, fbops, ffops = adaptDD(chan_op, ff_filt ,fb_filt, 1)
#quantizer_op, SD2 = adaptDD(chan_op, rls_ff ,rls_fb, 1)
bits = Int64[]
[push!(bits, (sign(real(training_seq[i]))==sign(real(quantizer_op[i]))) ? 1 : 0)  for i ∈ 1:min(length(training_seq),length(quantizer_op))] 
ber1= 1 - (sum(bits)/length(bits))
println(ber1)
|



# quantizer_op, error = adaptRLSTrain(chan_op, training_seq, ff_filt_len, fb_filt_len, 0.999)
# bits = Int64[]
# [push!(bits, (sign(real(training_seq[i]))==sign(real(quantizer_op[i]))) & (sign(imag(training_seq[i]))==sign(imag(quantizer_op[i]))) ? 1 : 0)  for i ∈ 1:min(length(training_seq),length(quantizer_op))] 
# ber1= 1 - (sum(bits)/length(bits))
# println(ber1)




