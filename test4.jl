include("rls.jl")

L = 6000
yb1 = yb[1:L]
D = 0

r=samples((bbc[events[1]:events[1]+L+D-1]))
r = resample(r,1000/999)
r=(r)/norm(r,2)*norm(yb1)
#r = r .+ 1
Q = 3
#r = r *norm(yb1,2) /norm(r,2)
quantizer_op, error,ff_filt, fb_filt, SD, ff1, ff2, aaa = adaptRLSTrain2(r, samples(yb1), 12*Q, 12*Q, 0.99, 1)
plot(real.(samples(yb[1:L])))
plot!(real.(SD))
bits = Int64[]
[push!(bits, (sign(real(yb[i]))==sign(real(SD[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(SD))]
ber1= 1 - (sum(bits)/length(bits))
println(ber1)
plot(angle.(samples(yb[1:L])))
plot!(abs.(angle.(SD)))
bits = Int64[]
atx=angle.(samples(yb[1:L]))
tx_m=mean(atx)
tx_zm = atx .- tx_m
aSD = abs.(angle.(SD))
aSD_zm = aSD .- tx_m
[push!(bits, (sign((tx_zm[i]))==sign((aSD_zm[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(SD))]
ber1= 1 - (sum(bits)/length(bits))
println(ber1)
quantizer_op, SD2, fbops, ffops = adaptDD(r, ff_filt, fb_filt, 1)
bits = Int64[]
aSD = abs.(angle.(SD2))
aSD_zm = aSD .- tx_m
[push!(bits, (sign((tx_zm[i]))==sign((aSD_zm[i])))  ? 1 : 0)  for i ∈ 1:min(length(yb),length(SD))]
ber1= 1 - (sum(bits)/length(bits))
println((ber1))

