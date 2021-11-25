
# y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
# yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)


# for i = 0:12
# rbb = bb[ events[1]+i-12000 : events[1] + T + i+ 12000]
# corr = xcorr(samples(rbb), samples(yb))
# #rbd = downsymbol(rbb)
# #corr = xcorr(rbd, mseq(12))
# display(plot(abs.(corr)))
# end

function downsymbol(yb)
  temp = ComplexF64[]
  for i = 12:12:length(yb)
    push!(temp,convert.(ComplexF64,mean(yb[i-11:i])))
  end
  temp
end



# for i = 1020:-1:950
#   if(i != 1000)
#     ss= sresample(rbb, 1000//i)
#     #sss= ss[1:length(ss)]
#     quantizer_op, error, sd1, ff_filt, fb_filt = adaptRLSTrain(samples(ss), txbb, 36, 12, 0.99, 1)
#   end
#   println((i,norm(error,2)))
# end


# rbb =  bb[ events[1]+6 : events[1] + T + 6] 
# ss= sresample(rbb, 1000//999)
# ss = ss*norm(txbb,2)/norm(ss,2)
# sss= ss[6:length(ss)]
# quantizer_op, error, sd1, ff_filt, fb_filt = adaptRLSTrain(samples(sss), txbb, 1, 1, 0.995, 1)
# println(norm(error,2)/length(error))

# for i = -3:0.1:3
#   for j = -3:0.1:3
#   bits = Int64[]
#   quantizer_op,  sd1, a1, a2 = adaptDD(samples(ss), ff_filt*j, fb_filt*i, 1)
#   [push!(bits, sign(real(txbb[i]))==sign(real(quantizer_op[i])) ? 1 : 0)  for i ∈ 1:min(length(txbb),length(quantizer_op))] 
#   ber1= 1 - (sum(bits)/length(bits))
#   println((i,j,ber1))  
#   end
# end



