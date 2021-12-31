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
    for j = 1:1:100
      vn2 = resample(vn1, (j+950)/1000)
      cr=mfilter(yb, vn2)
      max = maximum(abs.(cr))
      push!(a, max)
    end
  fac=findmax(a)[2]
  vn_best = resample(vn1, (fac+950)/1000)
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

function find_mseq(bb, tx, th, blk)
  λ = 0.9999          # exponential averaging factor for threshold
  β = 5.0             # threshold is β × average
  xmin = th       # minimum threshold
  gap = 3500         # minimum gap between detections
  pwidth = 25        # peak width (to look for maxima)
  x = abs.(mfilter(tx[1:blk], bb))
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
  #plot(x)
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
  if(length(invalid)!=1)
    deleteat!(events, findall(invalid))
  end
  events
end 

function process(bbc, events, d, Q)
  y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
  yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
  T = length(yb)
  rbb = bbc[events:events+T]
  vn = makevb(rbb, d) 
  vn1 = vn #resample(vn[1:min(length(yb), length(vn))],(1000)/1000)
  blksize = 1092
  eqlzBlkSize = (12 * 64)
  e2 = Float64[]
  vn_best = zeros(ComplexF64,1,blksize)
  aSD_zm = zeros(ComplexF64,1,blksize)
  tx_zm = zeros(ComplexF64,1,blksize)
  push!(e2 , 0.0)
  #length(vn1)
  λ = 0.9995 
  vn_best_all = ComplexF64[]
  SD_all = ComplexF64[]
  for i = 1:blksize:length(vn1)-blksize-1
      vn_best = compDoppler(vn1[i:i+blksize-1], yb[i:i+blksize-1])
      if(length(vn_best) > blksize)
          vn_best = vn_best[1:blksize]
      end
      if(length(vn_best) < blksize)
          vn_best = [vn_best; zeros(ComplexF64,blksize-length(vn_best))]
      end

      for h = 1:length(vn_best)
          push!(vn_best_all, vn_best[h])
      end
  end

  for i = 1:eqlzBlkSize:length(vn_best_all)-eqlzBlkSize-1
      quantizer_op, error, ff_filt, SD, ffops, ff = chRlsTrain(vn_best_all[i:i+eqlzBlkSize], yb[i:i+eqlzBlkSize], Q*12, λ, 1)
      Ber1 = computeBER(quantizer_op, yb[i:i+eqlzBlkSize])
      quantizer_op, SD2 = DA(vn_best_all[i:i+eqlzBlkSize], ff_filt, 1)
      println((length(SD2), eqlzBlkSize))
      for h = 1:(eqlzBlkSize)
        if(h <= length(SD2))
          push!(SD_all, SD2[h])
        else
          push!(SD_all, zero(ComplexF64))
        end
      end
  end



  atx=angle.(samples(yb))
  tx_m = mean(atx)
  aSD = abs.(angle.(SD_all))
  tx_zm = atx .- tx_m
  aSD_zm = aSD .- tx_m
  #println(length(aSD_zm))
  downbb, downtxbb = downsymboling(aSD_zm, tx_zm)
  ber = computeBER(downbb, downtxbb)

   
  # for i = 1:eqlzBlkSize:length(aSD_zm)

  # end


  println((length(aSD_zm),length(tx_zm),length(vn_best_all),ber))
  ber, downbb, downtxbb
end

function process3(events, d, Q, x)
  y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
  yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
  T = length(yb)
  # filename = "loc5-1.wav" 
  # nsamples, _ = wavsize(filename)
  # x = signal(filename; start = 1, nsamples = nsamples)
  bb1,_ = prepsig(x[:,1])
  bb2,_ = prepsig(x[:,2])
  bb3,_ = prepsig(x[:,3])
  println(events, length(bb1))
  rbb1 = bb1[events:events+T]
  rbb2 = bb2[events:events+T]
  rbb3 = bb3[events:events+T]

  vn1 = makevb(rbb1, d) 
  vn2 = makevb(rbb2, d)
  vn3 = makevb(rbb3, d)

  blksize = 300
  e2 = Float64[]
  vn_best = zeros(ComplexF64,1,blksize)
  aSD_zm = zeros(ComplexF64,1,blksize)
  tx_zm = zeros(ComplexF64,1,blksize)
  push!(e2 , 0.0)
  #length(vn1)
  λ = 0.9995 
  for i = 1:blksize:length(vn1)-blksize
      vn_best1 = compDoppler(vn1[i:i+blksize], yb[i:i+blksize])
      vn_best2 = compDoppler(vn2[i:i+blksize], yb[i:i+blksize])
      vn_best3 = compDoppler(vn3[i:i+blksize], yb[i:i+blksize])        
      if(length(vn_best2) > blksize)
          vn_best2 = vn_best2[1:blksize]
      end
      if(length(vn_best2) < blksize)
          vn_best2 = [vn_best2; zeros(ComplexF64,blksize-length(vn_best2))]
      end

      if(length(vn_best1) > blksize)
          vn_best1 = vn_best1[1:blksize]
      end
      if(length(vn_best1) < blksize)
          vn_best1 = [vn_best1; zeros(ComplexF64,blksize-length(vn_best1))]
      end

      if(length(vn_best3) > blksize)
          vn_best3 = vn_best3[1:blksize]
      end
      if(length(vn_best3) < blksize)
          vn_best3 = [vn_best3; zeros(ComplexF64,blksize-length(vn_best3))]
      end

      quantizer_op1, error, ff_filt1, SD, ffops, ff = chRlsTrain(vn_best1[1:blksize], yb[i:i+blksize], Q*12, λ, 1)
      Ber1 = computeBER(quantizer_op1, yb[i:i+blksize])
      quantizer_op2, error, ff_filt2, SD, ffops, ff = chRlsTrain(vn_best2[1:blksize], yb[i:i+blksize], Q*12, λ, 1)
      Ber2 = computeBER(quantizer_op2, yb[i:i+blksize])
      quantizer_op3, error, ff_filt3, SD, ffops, ff = chRlsTrain(vn_best3[1:blksize], yb[i:i+blksize], Q*12, λ, 1)
      Ber3 = computeBER(quantizer_op3, yb[i:i+blksize])
      println((i,Ber1, Ber2, Ber3))

      quantizer_op1, SD1 = DA(vn_best1[1:blksize], ff_filt1, 1)
      quantizer_op2, SD2 = DA(vn_best2[1:blksize], ff_filt2, 1)
      quantizer_op3, SD3 = DA(vn_best3[1:blksize], ff_filt3, 1)
      SDT = (0.3*SD1 + 0.5*SD2 + 0.2*SD3)
      atx=angle.(samples(yb[i:i+blksize]))
      tx_m=mean(atx)
      aSD = abs.(angle.(SDT))
      tx_zm = atx .- tx_m
      aSD_zm = aSD .- tx_m
      #println(length(aSD_zm))
      downbb, downtxbb = downsymboling(aSD_zm, tx_zm)
      ber = computeBER(downbb, downtxbb)
      println((i,ber))
      push!(e2,ber)
  end
  e2
end