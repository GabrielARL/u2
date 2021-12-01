using WAV
using SignalAnalysis
using DSP


const audacityExe = "/Applications/Audacity.app/Contents/MacOS/Audacity"
fs = 32000
bpf = fir(256, 3800.0, 6200.0; fs=fs)
bpf1 = fir(256, 4800.0, 5200.0; fs=fs)

function computeBER(downbb, downtxbb)
  bits = Int64[]
  [push!(bits, (sign((downbb[i]))==sign((downtxbb[i])))  ? 1 : 0)  for i ∈ 1:min(length(downbb),length(downtxbb))]
  ber2= 1 - (sum(bits)/length(bits))
  ber2
end

function downsymboling(aSD_zm, tx_zm)
  SDs = Float64[]
  atxDs = Float64[]
  downbb = Int64[]
  downtxbb = Int64[]


  for h = 1:length(aSD_zm)
    push!(SDs, aSD_zm[h])
    push!(atxDs, tx_zm[h])
    if(h%12==0)
      temp1 = sign.(aSD_zm[h:-1:(h-12+1)])
      n=count(>=(0.0), temp1)
      if (n >= 6)
        bit = 1
      else 
        bit = -1
      end
      push!(downbb,bit)
      temp1 = sign.(tx_zm[h:-1:(h-12+1)])
      n=count(>=(0.0), temp1)
      #println(n)
      if (n >= 6)
        bit = 1
      else 
        bit = -1
      end
      push!(downtxbb,bit)
    end
    
  end

  downbb, downtxbb

end


function compDoppler(vn1, yb)
    a = Float64[]
    for j = 1:1:200
      vn2 = resample(vn1, (j+900)/1000)
      cr=mfilter(yb, vn2)
      max = maximum(abs.(cr))
      push!(a, max)
    end
  fac=findmax(a)[2]
  vn_best = resample(vn1, (fac+900)/1000)
  vn_best
end

function prepsig(x)
    d1 = sfiltfilt(bpf, x)
    bb=downconvert(sresample(d1, 9//8), 6, 6000)
    bbc = bb .* cw(1000.0, length(bb)/6000, 6000.0)

    d1 = sfiltfilt(bpf1, x)
    bb=downconvert(sresample(d1, 9//8), 6, 6000)
    bbc2 = bb .* cw(1000.0, length(bb)/6000, 6000.0)
    bbc, bbc2
end


function audacity(x; fs=framerate(x), timeout=10.0)
  wavfname = tempname() * ".wav"
  wavwrite(samples(x), wavfname; Fs=fs)
  run(`"$audacityExe" "$wavfname"`; wait=false)
  @async begin
    sleep(timeout)
    rm(wavfname)
  end
  nothing
end

wavsize(filename) = wavread(filename; format="size")

function find_mseq(bb, tx)
  λ = 0.9999          # exponential averaging factor for threshold
  β = 5.0             # threshold is β × average
  xmin = 0.7       # minimum threshold
  gap = 3500         # minimum gap between detections
  pwidth = 25        # peak width (to look for maxima)
  x = abs.(mfilter(tx[1:3000], bb))
  μ = x[1]
  j = 0
  events = Int[]
  fs = 6000
  for i = 1:length(x)
    μ = λ * μ + (1 - λ) * x[i]
    if x[i] ≥ max(xmin, β * μ) && i ≥ j
      push!(events, argmax(x[i:min(length(x),i+pwidth)]) + i - 1)
      j = i + gap
    end
  end
  plot(x)
  println(events)

  invalid = ones(Bool, size(events))
  for i = 2:length(events)
    dt = (events[i] - events[i-1]) / fs
    n = dt/19  
    println(n)                                     # 19s gap between LFMs in a group
    if n > 1.0 && n < 3.0 && isapprox(n, round(n); atol=0.4)  # valid -> gap is right
      invalid[i] = false
      invalid[i-1] = false
    end
  end
  println(invalid)
  deleteat!(events, findall(invalid))
  events
end 