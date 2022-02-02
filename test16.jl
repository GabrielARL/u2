
# events = Int64[]
# events = find_mseq(bbc, (yb), 0.8, 3000) #0.7, 3000 for 4-2 #0.5, 3000 for 5-1 # 0.7, 2000 for 5-2 # 0.7,2000
#  pushfirst!(events, 228395)
#  pushfirst!(events, 98365)
#events = [105659, 224954, 371316, 471439, 563463+30000-989, 718312-177, 845395-4000+83, 966890+31000]#events 5-1

include("utils.jl")

filename = "loc6-1.wav" 
nsamples, _ = wavsize(filename)
x = signal(filename; start = 1, nsamples = nsamples)
yb = signal(repeat(mseq(12); inner=12), 6000)
T = length(yb)
bb1, _ = prepsig(x[:,1])
bb2, _ = prepsig(x[:,2])
bb3, _ = prepsig(x[:,3])


events = find_mseq(bbc, (yb), 0.5, 3000)
events_calibed = Int64[]

for k = 1:length(events)
a = Int64[]
b = Float64[]
for i = -30000:1000:30000
snr1,snr2,snr3 = get_snr(events[k]+i, bb1, bb2, bb3, T)
#println((snr1,snr2,snr3,i))
push!(a, i)
push!(b, snr2)
end
(v,q)=findmax(b) 
maxx=a[q]
a = Int64[]
b = Float64[]
for i = -1000:1:1000
snr1,snr2,snr3 = get_snr(events[k]+maxx+i, bb1, bb2, bb3, T)
#println((snr1,snr2,snr3,i))
push!(a, i)
push!(b, snr2)
end
(v,q)=findmax(b) 
max1=a[q]
push!(events_calibed, maxx+max1)
end

events = events .+ events_calibed