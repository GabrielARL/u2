using WAV
using SignalAnalysis
using DSP


const audacityExe = "/Applications/Audacity.app/Contents/MacOS/Audacity"
fs = 32000
bpf = fir(256, 3800.0, 6200.0; fs=fs)
bpf1 = fir(256, 4500.0, 5500.0; fs=fs)

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

function dSample_comp_err(SD_all, yb)
  atx=angle.(samples(yb))
  tx_m = mean(atx)
  aSD = abs.(angle.(SD_all))
  tx_zm = atx .- tx_m
  aSD_zm = aSD .- tx_m
  #println(length(aSD_zm))
  downbb, downtxbb = downsymboling(aSD_zm, tx_zm)
  uber = computeBER(downbb, downtxbb)
  OSNR=10*log10(1/(sum((abs.((yb[1:length(SD_all)] .- (-SD_all) ))).^2)/length(SD_all)))
  uber, downbb, downtxbb, OSNR
end

function switchBool(a)
  if(a == -1)
    b = 0
  else
    b = 1
  end
  convert(Bool, b)
  b
end

function perform_post_decoding(data, downbb,downtxbb)
    EncData=FEC.encode(BCH.BCH_31_6, data)
    sc = xor.(EncData, downtxbb[1:length(EncData)])

    EncDataR1 = xor.(sc, downbb[1:length(EncData)])
    EncDataR_buffer = ones(Int64,1,length(EncDataR1))
    for i = 1:length(EncDataR1)
        if(EncDataR1[i] == 0)
        EncDataR_buffer[i] = -1
        end
    end
  dataR , _ = FEC.decode(BCH.BCH_31_6, EncDataR_buffer)
  dataR
end



function process3(x, events, d, Q)
  y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
  yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)
  T = length(yb)

  bb1, _ = prepsig(x[:,1])
  bb2, _ = prepsig(x[:,2])
  bb3, _ = prepsig(x[:,3])

  rbb1 = bb1[events:events+T]
  rbb2 = bb2[events:events+T]
  rbb3 = bb3[events:events+T]

  nbb1 = bb1[events+T:min(events+(2*T),length(bb1)-1)]  
  nbb2 = bb2[events+T:min(events+(2*T),length(bb1)-1)]
  nbb3 = bb3[events+T:min(events+(2*T),length(bb1)-1)]

  snr1 = norm(rbb1,2) / ( norm(nbb1, 2))
  snr2 = norm(rbb2,2) / ( norm(nbb2, 2))
  snr3 = norm(rbb3,2) / ( norm(nbb3, 2))

  println((snr1,snr2,snr3))

  vn1 = makevb(rbb1, d)
  vn2 = makevb(rbb2, d)
  vn3 = makevb(rbb3, d) 
  
  blksize = 1092
  eqlzBlkSize = (12 * 64)
  e2 = Float64[]
  push!(e2 , 0.0)
  #length(vn1)
  λ = 0.9995 

  SD_all1 = ComplexF64[]
  SD_all2 = ComplexF64[]
  SD_all3 = ComplexF64[]

  data = BitVector(rand(Bool, 792))
  SD_all1 = compDoppler_equalize(vn1, blksize, eqlzBlkSize, yb, Q, λ)
  SD_all2 = compDoppler_equalize(vn2, blksize, eqlzBlkSize, yb, Q, λ)
  SD_all3 = compDoppler_equalize(vn3, blksize, eqlzBlkSize, yb, Q, λ)

  uber1, downbb1, downtxbb, OSNR1 = dSample_comp_err(SD_all1, yb)
  uber2, downbb2, downtxbb, OSNR2 = dSample_comp_err(SD_all2, yb)
  uber3, downbb3, downtxbb, OSNR3 = dSample_comp_err(SD_all3, yb)
  
  OSNRT = OSNR1 + OSNR2 + OSNR3
  SD_all_c = ((SD_all1 * OSNR1) + (SD_all2*OSNR2) + (SD_all3*OSNR3))/  OSNRT
  uber4, downbb4, downtxbb, OSNR4 = dSample_comp_err(SD_all_c, yb)

  downbb1 = [downbb1; ones(ComplexF64,length(mseq(12))-length(downbb1))]
  downbb2 = [downbb2; ones(ComplexF64,length(mseq(12))-length(downbb2))]
  downbb3 = [downbb3; ones(ComplexF64,length(mseq(12))-length(downbb3))]
  downbb4 = [downbb4; ones(ComplexF64,length(mseq(12))-length(downbb4))]
  downtxbb = [downtxbb; ones(ComplexF64,length(mseq(12))-length(downtxbb))]
  
   downbb1=switchBool.(downbb1)
   downbb2=switchBool.(downbb2)
   downbb3=switchBool.(downbb3)
   downbb4=switchBool.(downbb4)
   downtxbb = switchBool.(downtxbb)

   dataR1 = perform_post_decoding(data, downbb1, downtxbb) 
   dataR2 = perform_post_decoding(data, downbb2, downtxbb)
   dataR3 = perform_post_decoding(data, downbb3, downtxbb)
   dataR4 = perform_post_decoding(data, downbb4, downtxbb)

  cber1 = sum(abs.(data .- dataR1))/length(data)
  cber2 = sum(abs.(data .- dataR2))/length(data) 
  cber3 = sum(abs.(data .- dataR3))/length(data)
  cber4 = sum(abs.(data .- dataR4))/length(data)
  cber1, cber2, cber3, cber4
end