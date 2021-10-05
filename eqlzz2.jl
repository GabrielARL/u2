using LinearAlgebra
using DSP
Random.seed!(1234);

include("gabUtil.jl")

function chsign(a)
    if(a == 0)
        return -1
    end
    return 1
end

training_len = 10^4
snr_dB = 30
ff_filt_len = 2
fb_filt_len = 2
snr = 10^(0.1*snr_dB)
noise_var = 1/(2*snr)

training_a = rand([1 0], training_len)
training_seq = 1 .- 2 .* training_a .+0*im; 

fade_chan = [0.407+0.01im 0.815+0.02im 0.1+0.01im]'; # Proakis B channel
fade_chan = fade_chan/norm(fade_chan,2);
chan_len = length(fade_chan);
noise = sqrt(noise_var).*randn(1, training_len+chan_len-1)
chan_op = conv(fade_chan, training_seq) .+ noise'

ff_filt = zeros(1,ff_filt_len); 
fb_filt = zeros(1,fb_filt_len); 
ff_filt_ip = zeros(1,ff_filt_len); 
fb_filt_ip = zeros(1,fb_filt_len); 

fb_filt_op = 0;
Rvv0 = (chan_op'*chan_op)/(training_len+chan_len-1);
max_step_size = 2 ./ (ff_filt_len*Rvv0.+fb_filt_len*(1));
step_size = 0.125*max_step_size; # step size

quantizer_op, ff_filt, fb_filt, err = adaptLMSTrain(chan_op, training_seq, ff_filt_len, fb_filt_len, 0.001)

bits = Int[]
dec_a = sign.(quantizer_op) 
dec = convert.(Int64, dec_a)
data = convert.(Int64, training_seq)
#dec_a = dec_a[1:training_len-ff_filt_len+1]
#data = chsign.(data)
[push!(bits, data[i]==dec[i] ? 1 : 0)  for i ∈ 1000:length(dec)] 
ber=1-(sum(bits)/length(bits))



# data_len = training_len
# data_a = rand([1 0], training_len)
# data_seq = 1 .- 2 .* data_a;
# noise = sqrt(noise_var).*randn(1, training_len+chan_len-1)
# chan_op = conv(fade_chan, data_seq) .+ noise'

# dec_seq = zeros(1,data_len-ff_filt_len+1);

# ff_filt_ip = zeros(1,ff_filt_len); 
# fb_filt_ip = zeros(1,fb_filt_len); 
# fb_filt_op = 0;

# for i1=1:data_len-ff_filt_len+1 
#     ff_filt_ip[2:end]=ff_filt_ip[1:end-1];
#     ff_filt_ip[1] = chan_op[i1];
#     ff_filt_op = vec(ff_filt)'*vec(ff_filt_ip)
#     ff_and_fb = ff_filt_op-fb_filt_op; 
#     temp = ff_and_fb<0;
#     dec_seq[i1] = 1-2*temp;

#     fb_filt_ip[2:end]=fb_filt[1:end-1];
#     fb_filt_ip[1] = dec_seq[i1];
#     fb_filt_op = vec(fb_filt)'*vec(fb_filt_ip);
# end

# bits = Int[]
# dec_a = dec_seq.<0
# dec = convert.(Int64, dec_a)
# data = convert.(Int64, data_a)
# dec_a = dec_a[1:training_len-ff_filt_len+1]
# [push!(bits, data[i]==dec[i] ? 1 : 0)  for i ∈ eachindex(dec)] 
# ber=1-(sum(bits)/length(bits))