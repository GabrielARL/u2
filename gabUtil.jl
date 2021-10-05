using Statistics
using DataFrames
using Query
using LinearAlgebra
using CSV
using SignalAnalysis
using DSP
using Plots
using FileIO
using StatsPlots
using Printf
using Dates
using Geodesy
using ToeplitzMatrices

const base = @__DIR__
const fs = 1.0 / (1.6e-6 * 26)
const rx2µPa = db2amp(148)
const good = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] 
const bpf = fir(256, 3800.0, 6200.0; fs=fs)
const calib = [1.0 1.0 1.0 1.0 0.94291 1.0 1.0 1.0 0.988606 1.0 1.04269 0.76304] 
const lfm = real(chirp(4000.0, 6000.0, 1.0, fs; window=(tukey, 0.05)))
const β = 7
const xmin = 1000
const Δ1 = 0.1

function loadsoundings()
    data = NTuple{3,Float64}[]
    for filename ∈ filter(s -> endswith(s, ".txt"), readdir(joinpath(base, "data/ENC"), join=true))
      open(filename, "r") do io
        for s ∈ eachline(io)
          m = match(r"^<SOUNDG .*>\( (.*) \)", s)
          if m !== nothing
            for x ∈ split(m[1], " ), ( ")
              try
                y = parse.(Float64, split(x, ","))
                push!(data, (y[1], y[2], y[3]))
              catch
                @warn "Bad data: $s"
              end
            end
          end
        end
      end
    end
    data
  end

  function allrecordings()
    out = DataFrame(id=Int64[], dts=DateTime[], distance=Float64[], srclvl=Float64[], txdepth=Float64[], rxdepth=Float64[], filename=String[], duration=Float64[])
    logs = CSV.read(joinpath(base, "data/recinfo.csv"))
    for r ∈ eachrow(logs)
      id, dd, hhmmss, rx2mins, tx2mins, rxlat1, rxlon1, rxlat2, rxlon2, txlat1, txlon1, txlat2, txlon2, txdepth, srclvl, anchor = r
      txlat1 == 0 && txlon1 == 0 && continue
      r = distance(LLA(rxlat1, rxlon1), LLA(txlat1, txlon1))
      dts = DateTime(@sprintf("202009%02d%06d", dd, hhmmss), "yyyymmddHHMMSS")
      filename = joinpath(base, @sprintf("data/202009%02d/09%02d20_%06d.dat", dd, dd, hhmmss))
      info = stat(filename)
      nsamples = floor(Int, info.size / 52)
      if nsamples > 0 && (srclvl > 170 || r < 10000.0)
        push!(out, (id, dts, r, srclvl, txdepth, anchor, filename, nsamples/fs))
      end
    end
    out[ 7,:duration] -= 185.0
    out[10,:duration] -=  55.0
    out[11,:duration] -= 230.0
    out[12,:duration] -= 185.0
    out[13,:duration]  = 150.0
    delete!(out, [8,9])
    out
  end

  function readdtla(filename, start=0.0, duration=nothing)
    framelen = 2 * 26
    channels = 24
    magic = 0xc0de
    info = stat(filename)
    info.size >= 4*channels || error("File too small")
    nsamples = floor(Int, info.size / framelen)
    t = nsamples / fs
    first = round(Int, start * fs)
    nsamples -= first
    nsamples <= 0 && error("Invalid requested start")
    if duration !== nothing
      n = round(Int, duration * fs)
      nsamples < n && (n = nsamples)
      nsamples = n
    elseif nsamples > 60 * fs
      @warn "Loading 1 minute of data (total $(round(Int, t/6)/10) minutes)..."
      nsamples = round(Int, 60 * fs)
    end
    open(filename, "r") do f
      data = Array{UInt16}(undef, framelen >> 1, nsamples)
      seek(f, round(Int, first * framelen))
      read!(f, data)
      data[1,1] == magic || error("Invalid file format")
      data[2,1] == magic || error("Invalid file format")
      x = @. 5.0 * (data[3:end,:]' / 65536.0) - 2.5
      signal(x * rx2µPa, fs)
    end
  end

  function lfm_amplitude(s, i)
    s1 = @view s[max(1, i - 250) : min(length(s), i + length(lfm) + 250)]
    θ = 4 * std(s1)
    maximum(abs.(mfilter(lfm, clamp.(s1, -θ, θ)))[1:501])
  end

  function find_lfms2(s)
    λ = 0.9999          # exponential averaging factor for threshold #mandar: 0.9999
    #println((β, xmin))
    β = 5
    xmin = 700
    gap = 25000         # minimum gap between detections
    pwidth = 250        # peak width (to look for maxima)
    x = abs.(mfilter(lfm, sign.(s)))
    μ = x[1]
    j = 0
    events = Int[]
    for i = 1:length(x)
      μ = λ * μ + (1 - λ) * x[i]
      if x[i] ≥ max(xmin, β * μ) && i ≥ j
        push!(events, argmax(x[i:min(length(x),i+pwidth)]) + i - 1)
        j = i + gap
      end
    end
    invalid = ones(Bool, size(events))
    for i = 2:length(events)
      dt = (events[i] - events[i-1]) / fs
      n = dt/2.6                                        # 2.6s gap between LFMs in a group
      if dt < 10.0 && isapprox(n, round(n); atol=0.01)  # valid -> gap is right
        invalid[i] = false
        invalid[i-1] = false
      end
    end
    deleteat!(events, findall(invalid))
    E = sum(lfm .^ 2)
    DataFrame(
      time = events,
      detector = samples(x[events] ./ sum(abs.(lfm))),
      rxlvl = [amp2db(lfm_amplitude(s, i)/E) for i ∈ events]
    )
  end

  function find_lfms(s)
    λ = 0.9999          # exponential averaging factor for threshold #mandar: 0.9999
    #println((β, xmin))
    #β = 7
    #xmin = 1000
    gap = 25000         # minimum gap between detections
    pwidth = 250        # peak width (to look for maxima)
    x = abs.(mfilter(lfm, sign.(s)))
    μ = x[1]
    j = 0
    events = Int[]
    for i = 1:length(x)
      μ = λ * μ + (1 - λ) * x[i]
      if x[i] ≥ max(xmin, β * μ) && i ≥ j
        push!(events, argmax(x[i:min(length(x),i+pwidth)]) + i - 1)
        j = i + gap
      end
    end
    invalid = ones(Bool, size(events))
    for i = 2:length(events)
      dt = (events[i] - events[i-1]) / fs
      n = dt/2.6                                        # 2.6s gap between LFMs in a group
      if dt < 10.0 && isapprox(n, round(n); atol=0.01)  # valid -> gap is right
        invalid[i] = false
        invalid[i-1] = false
      end
    end
    deleteat!(events, findall(invalid))
    E = sum(lfm .^ 2)
    DataFrame(
      time = events,
      detector = samples(x[events] ./ sum(abs.(lfm))),
      rxlvl = [amp2db(lfm_amplitude(s, i)/E) for i ∈ events]
    )
  end

  function process_mseq(ndx, groups, data)
    Δ = round(Int, (5*2.6+10)*fs)
    s = gmseq(8)
    T = round(Int, 5*length(s)*128*fs/44100+2fs)
    df = DataFrame(time = [], rxlvl = [], group = [], count = [], doppler = [], ir = [])
    for n ∈ 1:length(ndx)
      i = ndx[n]
      group = groups[n]
      if i+Δ+T < length(data)
        x = sresample(data[i+Δ-1000:i+Δ+T], 44100/fs)
        y = downconvert(x, 128, 5000.0, rcosfir(0.25, 128))
        Δflist = -10:0.5:10
        z = Array{Float64}(undef, length(y), length(Δflist))
        for j ∈ 1:length(Δflist)
          Δf = Δflist[j]
          z[:,j] = abs.(mfilter(s, y .* cw(Δf, duration(y), framerate(y)))) / length(s)
        end
        z̄ = maximum(z; dims=2)[:,1]
        p = 6
        for j ∈ 1:5
          fp = firstpeak(z̄[p:min(p+1100,length(z̄))])
          #println(fp)
          fp === nothing && break
          p += fp - 1
          p1 = p + argmax(z̄[p:min(p+50,length(z̄))]) - 1
          k = argmax(sum(z[p1-1:p1+1,:]; dims=1)[:])
          z1 = z[p-5:min(p+length(s)-5,size(z,1)), k]
          push!(df, (i+Δ-1000 + round(Int, fs/44100*p), amp2db(z[p]), group, j, Δflist[k], z1))
          p += length(s) - 10
        end
      end
    end
    df
  end

  function firstpeak(x)
    #println(length(x))
    length(x) == 0 && return nothing
    θ = max(maximum(x)/3, 1.5*median(x))
    p = findfirst(x .≥ θ)
    #println(p)
    while p !== nothing && p < length(x) && x[p+1] > x[p]
     p += 1
     #println(p)
    end
    p
  end

  function label_lfms!(lfms)
    labels = Array{Int}(undef, size(lfms,1), 2)
    labels[1,:] .= [1, 1]
    for i ∈ 2:size(labels,1)
      dt = (lfms[i,:time] - lfms[i-1,:time]) / fs
      
      if ( dt > 2.5 && dt < 2.7)
        labels[i,1] = labels[i-1,1]         # 2.6s gap between LFMs in a group
        #println(dt)
        labels[i,2] = labels[i-1,2] + round(Int, dt/2.6)
      else
        n = round(Int, dt/65.0)             # 65s signal repeat 
        #println((dt, n, i, round(Int, (dt - 65.0 * n) / 2.6)))
        labels[i,1] = labels[i-1,1] + n
        labels[i,2] = labels[i-1,2] + round(Int, (dt - 65.0 * n) / 2.6)
      end
    end
    lfms[:group] = labels[:,1]
    lfms[:count] = labels[:,2] .- minimum(labels[:,2]) .+ 1
    #println(labels)
    lfms
  end

  function getTimings(recNo)
    recs = allrecordings()
    data = readdtla(recs[recNo,:].filename, 0, recs[recNo,:].duration)
    data = data[:,good]
    data = sfiltfilt(bpf, data)
    data .*= calib
    timings = DataFrame(start = [], ending = [], group = [], ch = [])
    for i ∈ 1:length(good)
      if(recNo != 8)
      lfms = find_lfms(@view data[:,i])
      else
        lfms = find_lfms2(@view data[:,i])
      end
      label_lfms!(lfms)
      lfms1 = lfms[lfms.count .== 1,:]
      mseq = process_mseq(lfms1.time, lfms1.group, @view data[:,i])
      if(size(mseq,1) >=1)
        for j = 1:maximum(mseq[:group])
          mseqfil = filter(row -> row.group == j, mseq)
          if(size(mseqfil, 1) != 0)
            maxcnt = maximum(mseqfil[:count])
            start = mseqfil[:time][1]
            ending = mseqfil[:time][maxcnt]+round(Int, fs/44100*length(gmseq(8)))
            push!(timings, (start, ending, j, i))
          end
        end
      end
    end
    return data, timings
  end

  #This is the MMSE Equalizer  
  function MMSE_eqlz(s, MSEQ, M)
    ϕff = xcorr(samples(s), samples(s))
    ϕfd = xcorr(MSEQ, samples(s))
    N = length(MSEQ)
    𝑅ff = ϕff[N:N+M-1]
    𝑅 = Toeplitz(𝑅ff, vec(𝑅ff))
    P = ϕfd[N:N+M-1]
    h = inv(collect(𝑅)) * (P)
    h
  end

#pack the signals into a data structure with 3 dimensions. first dimension is the group which is the the whole packet train sequence
#second is the number of elements, third is the samples itself
function  packSignals(data, timings)
  dataDf = Array{Union{Float64, Nothing}}(nothing,maximum(timings[:group]),length(good),90495)
  for i = 1:maximum(timings[:group])
    dff = filter(row->row.group==i, timings)
    if(size(dff,1) != 0)
      for j = 1:maximum(dff[:ch])
        dfff = filter(row->row.ch==j && row.group ==i, dff)
        if(size(dfff,1) != 0)
          diff = dfff[:ending][1] - dfff[:start][1]
          if( diff >= 50 && diff <= 1200 )
            dd = data[dfff[:start][1]:dfff[:start][1]+90494,j]
            mseqsig = dd/maximum(dd)
            dataDf[i,j,:] = mseqsig
          end
        end
      end
      
    end
  end
  dataDf, timings
end

#this is a wrapper to convert the signal and pass it to the Equalizer trainer 
function eqlzz(sig)
  if(!isnothing(sig[1]))
    sig1=convert.(Float64, sig)
    bb = repeat(gmseq(8), 5)
    x = sresample(sig1, 44100/fs)
    x = x[1:end-1]
    N_ff = 10
    N_fb = 2
    λ = 0.9
    eqlz,w,wb,e=adaptRLSTrain(x, bb, N_ff, N_fb, λ, 0.01)
    println(abs.(eqlz))
    #softDecision = adaptRLSTrain(samples(x),bb, N_ff, N_fb, λ, 0.01)
  end
  SoftDecision=eqlz./(mean(abs.(eqlz)))#*cis(-5/180*2π)
  SoftDecision
end

#takes in the data, returns the SNR for all group of that channel - return 2d array
function combineSoftDecisions1channel(dataDf, timings)
  SDs = Array{Union{ComplexF64,Nothing}}(nothing, (maximum(timings[:group]), length(good), 1272))
  OSNR = Array{Union{Float64,Nothing}}(nothing, (maximum(timings[:group]), length(good)))
  for i = 1:maximum(timings[:group])
    dff = filter(row->row.group==i, timings)
    if(size(dff,1) != 0)
        for j ∈ good
        if(!isnothing(dataDf[i,j,1]))
          SDs[i,j,:] = eqlzz(dataDf[i, j,:])
          OSNR[i,j]=compOSNR(SDs[i,j,:])
        end
      end
    end
  end
  OSNR
end

#to compute OSNR
function compOSNR(s)
  bb = repeat(gmseq(8), 5)
  OSNR = 0
  if(!isnothing(s[1]))
  OSNR=10*log10(1/(sum((abs.((bb[80:length(s)] .- (-s[80:end]) ))).^2)/length(s[80:end])))
  end
  OSNR
end

#EGC combining of n channels,returns OSNRs, 
#the elements in each ith row and jth column 
#represents the OSNR of the i,j 2 channel combination pair 
function combineSoftDecisionsAllChannels(dataDf, timings)
  OSNRs = Array{Union{Float64,Nothing}}(nothing, (maximum(timings[:group]), length(good), length(good)))
  for i = 1:maximum(timings[:group])
    for j = 1:length(good)
      for k = 1:length(good)
          if(j<k)
            if(!isnothing(dataDf[i,j,1]) && !isnothing(dataDf[i,k,1]))
              OSNRs[i,j,k]  =  compOSNR((eqlzz(dataDf[i,j,:]) + eqlzz(dataDf[i,k,:]) ) /2)
            end
          end
      end
    end
  end           
  OSNRs     
end

#computes the avegare of OSNR across each channel for all group. returns a vector with size of length of channels. 
function aveOSNR(OSNRs)
  aveOSNRs=Array{Union{Float64,Nothing}}(nothing, length(good), length(good))
  aveVars=Array{Union{Float64,Nothing}}(nothing, length(good), length(good))
  for i = 1:length(good)
    for j = 1:length(good)
       mn, vars= mean_var_CFNothing(OSNRs[:,i,j])
       aveOSNRs[i,j] = mn[1]
       aveVars[i,j] = vars[1]
    end
  end
  aveOSNRs, aveVars
end


#collate all the OSNRs according to thier separation distance and compute their mean and variance
function aveOSNR_dist(aveOSNRs, aveVars)
  OSNR_Dist=Array{Float64}(undef, length(good)-1)
  OSNRvar_Dist=Array{Float64}(undef, length(good)-1)

  bb = Float64[]

  for i = 2:length(good)
  aa=  (triu(aveOSNRs,i-1)-triu(aveOSNRs,i))
     empty!(bb)
     for k ∈ CartesianIndices(aa)
          if(!isnothing(aa[k]) &&  aa[k] > 0.0 && !isnan(aa[k]))
            push!(bb, aa[k])
          end 
     end
     OSNR_Dist[i-1] =  -sum(-1 .* bb)/count(!iszero,bb)
  end

  for i = 2:length(good)
    aa = (triu(aveVars,i-1)-triu(aveVars,i))
    empty!(bb)
    for k ∈ CartesianIndices(aa)
      if(!isnothing(aa[k]) &&  aa[k] > 0.0 && !isnan(aa[k]))
        push!(bb, aa[k])
      end 
 end
    OSNRvar_Dist[i-1] =  -sum(-1 .* bb)/count(!iszero,bb)
  end
  OSNR_Dist, OSNRvar_Dist
end
  
function computeOSNRonly(combined_results)
  bb = repeat(gmseq(8), 5)
  OSNR = Array{Union{Float64,Nothing}}(nothing,size(combined_results,1), size(combined_results,2))
  for i = 1:size(combined_results,1)
    for j = 1:size(combined_results,2)
      if(!isnothing(combined_results[i,j,1]))
        OSNR[i,j]=10*log10(1/(sum((abs.((bb[1:length(combined_results[i,j,:])] - (-combined_results[i,j,:]) ))).^2)/length(combined_results[i,j,:])))
      else
        #println((i,j))  
      end
    end
  end

  contain_osnr = Float64[]
  ave_osnr = Array{Float64}(undef,size(OSNR,2)) 
  var_osnr = Array{Float64}(undef,size(OSNR,2)) 
  for i = 1:size(OSNR,2)
  cnt = 0
  for j = 1:size(OSNR,1)
    if(!isnothing(OSNR[j,i]))
      push!(contain_osnr, convert.(Float64,OSNR[j,i]))
      ave_osnr[i]+=OSNR[j,i]
      cnt +=1
      #println((j,i, OSNR[j,i]))
    end
  end
  if(cnt!=0)
    ave_osnr[i] = ave_osnr[i]/cnt
  end

  var_osnr[i] = var(contain_osnr)
  end
  OSNR, ave_osnr
end

# to get the mean and var for a type of an union if Float64 and nothing
FloatOrNothingArray = Array{Union{Float64,Nothing}}
function mean_var_CFNothing(a::FloatOrNothingArray)
    means = Float64[]
    vars = Float64[]
    container = Float64[]

    for i ∈ 1:size(a,2)
      for j ∈ 1:size(a,1)
          if(!isnothing(a[j,i]))
            temp = convert(Float64, a[j,i])
            push!(container, temp) 
          end
      end
      mn = mean(container)
      var1 = var(container)
      empty!(container)
      push!(means, mn)
      push!(vars, var1)
    end
    means, vars
end

function computeOSNR(combined_results)
  bb = repeat(gmseq(8), 5)
  OSNR = Array{Union{Float64,Nothing}}(nothing,size(combined_results,1), size(combined_results,2))
  for i = 1:size(combined_results,1)
    for j = 1:size(combined_results,2)
      if(!isnothing(combined_results[i,j,1]))
        OSNR[i,j]=10*log10(1/(sum((abs.((bb[1:length(combined_results[i,j,:])] - (-combined_results[i,j,:]) ))).^2)/length(combined_results[i,j,:])))
      else
        #println((i,j))  
      end
    end
  end

  contain_osnr = Float64[]
  ave_osnr = Array{Float64}(undef,size(OSNR,2)) 
  var_osnr = Array{Float64}(undef,size(OSNR,2)) 
  for i = 1:size(OSNR,2)
  cnt = 0
  for j = 1:size(OSNR,1)
    if(!iszero(OSNR[j,i]))
      push!(contain_osnr, convert.(Float64,OSNR[j,i]))
      ave_osnr[i]+=OSNR[j,i]
      cnt +=1
      #println((j,i, OSNR[j,i]))
    end
  end
  if(cnt!=0)
    ave_osnr[i] = ave_osnr[i]/cnt
  end

  var_osnr[i] = var(contain_osnr)
  end
  OSNR,ave_osnr,var_osnr
end

  #this is the adaptive LMS DFE  
  function adaptLMSTrain(s, bb, N_ff, N_fb, δ)
    ff_filt_len = N_ff
    fb_filt_len = N_fb
    ff_filt = zeros(ComplexF64,ff_filt_len)
    fb_filt = zeros(ComplexF64,fb_filt_len)
    ff_filt_ip = zeros(ComplexF64,ff_filt_len)
    fb_filt_ip = zeros(ComplexF64,fb_filt_len)
  
    fb_filt_op = 0
    quantizer_op = zeros(ComplexF64,(length(bb)-ff_filt_len+1))
    err = zeros(ComplexF64,(length(bb)-ff_filt_len+1))
    for i1=1:length(bb)-ff_filt_len+1 
      ff_filt_ip[2:end]=ff_filt_ip[1:end-1];
      ff_filt_ip[1] = s[i1];
      ff_filt_op = (ff_filt)'*(ff_filt_ip)
      ff_and_fb = ff_filt_op-fb_filt_op; 
      err[i1] = (ff_and_fb - bb[i1])  
      
      #println((i1, fb_filt_op))
      #temp = imag(ff_and_fb) < 0.0 ? 0 : 1
      quantizer_op[i1] = sign(real(ff_and_fb)) + 0im
      ff_filt=ff_filt-δ*err[i1]*conj.(ff_filt_ip);
      fb_filt=fb_filt+δ*err[i1]*conj.(fb_filt_ip);

      if(any(i -> i > 10, abs.(ff_filt)) || abs(err[i1]) > 10)
        println((i1 ,(ff_filt_op), (err[i1]), (ff_and_fb)))
        println(fb_filt, fb_filt_ip)
        println(fb_filt_op)
        break
      end

      fb_filt_ip[2:end]=fb_filt[1:end-1];
      fb_filt_ip[1] = quantizer_op[i1];
      fb_filt_op = (fb_filt)'*(fb_filt_ip);
      if(any(isinf,fb_filt_op))
        println((i1, err[i1], fb_filt_op, any(isnan,fb_filt), any(isnan,ff_filt_ip)))
        break
      end 
    end
    collect(vec(quantizer_op)), ff_filt, fb_filt, collect(vec(err))
  end

  function adaptRLSTrain(s, bb, N_ff, N_fb, λ, Δ)
    ff_filt_len = N_ff
    fb_filt_len = N_fb
  
    ff_filt = zeros(ComplexF64,ff_filt_len)
    fb_filt = zeros(ComplexF64,fb_filt_len)
    ff_filt_ip = zeros(ComplexF64,ff_filt_len)
    fb_filt_ip = zeros(ComplexF64,fb_filt_len)
  
    fb_filt_op = 0
    quantizer_op = zeros(ComplexF64, length(s)-ff_filt_len-fb_filt_len)
    SD = zeros(ComplexF64, length(s)-ff_filt_len-fb_filt_len)
    P1 = Δ * I(N_ff)
    P2 = Δ * I(N_fb)

    err = zeros(ComplexF64, length(s)-ff_filt_len-fb_filt_len)
    #println(length(s)-ff_filt_len-fb_filt_len)
    for i1=1:length(s)-ff_filt_len-fb_filt_len
      ff_filt_ip[2:end]=ff_filt_ip[1:end-1]
      ff_filt_ip[1] = s[i1]
      ff_filt_op = (ff_filt)'*(ff_filt_ip)
      ff_and_fb = ff_filt_op-fb_filt_op; 
      err[i1] =  -(bb[i1] - ff_and_fb) 
      #println((i1, fb_filt_op))
      #temp = imag(ff_and_fb) < 0.0 ? 0 : 1
      quantizer_op[i1] = sign(real(ff_and_fb))
      SD[i1] = ff_and_fb
      #println((i1, ff_and_fb))
      u1 = ff_filt_ip
      nume1 = P1*u1
      denom1 = λ + u1'*nume1
      K1 = nume1/denom1
      PPrime1 = K1*u1'*P1
      PPrime1 = collect(PPrime1)
      P1=(P1.-PPrime1)/λ
      
      u2 = vec(fb_filt_ip)
      nume2 = P2*u2
      denom2 = λ + u2'*nume2
      K2 = nume2/denom2
      PPrime2 = K2*u2'*P2
      PPrime2 = collect(PPrime2)
      if( N_fb == 1 )
        P2= (only(P2)-only(PPrime2)) /λ
      else
        P2=(P2-(PPrime2))/λ
      end
      ff_filt=ff_filt-K1*err[i1];
      if( N_fb == 1)
        fb_filt=only(fb_filt)+only(K2*err[i1]);
        fb_filt_ip[1] = quantizer_op[i1];
        fb_filt_op = only(fb_filt * fb_filt_ip)
      else
        fb_filt=fb_filt+K2*err[i1];
        fb_filt_ip[2:end]=fb_filt[1:end-1];
        fb_filt_ip[1] = quantizer_op[i1];
        fb_filt_op = vec(fb_filt)'*vec(fb_filt_ip);
      end
      #fb_filt_op = vec(fb_filt)'*vec(fb_filt_ip);
    end
    quantizer_op[1:end], ff_filt, fb_filt, err, SD
  end

#This is the adaptive RLS Equalizer with feedback in decision directed mode (Todo: implement PLL)
function adaptRLSDD(s, ff_filt ,fb_filt)
  ff_filt_len = length(ff_filt)
  ff_filt_ip = zeros(1,length(ff_filt)); 
  fb_filt_ip = zeros(1,length(fb_filt)); 
  fb_filt_op = 0;
  quantizer_op = Array{ComplexF64}(undef,(length(s)-length(ff_filt)))
  SD = Array{ComplexF64}(undef,(length(s)-ff_filt_len))
  for i1=1:length(s)-length(ff_filt)
    ff_filt_ip[2:end]=ff_filt_ip[1:end-1];
    ff_filt_ip[1] = s[i1];
    ff_filt_op = vec(ff_filt)'*vec(ff_filt_ip)
    ff_and_fb = ff_filt_op-fb_filt_op; 
    SD[i1] = ff_and_fb
    quantizer_op[i1] = sign(real(ff_and_fb))
    if(length(fb_filt)==1)
    else
    fb_filt_ip[2:end]=fb_filt[1:end-1];
    end
    fb_filt_ip[1] = quantizer_op[i1];
    if( length(fb_filt)==1 )
      fb_filt_op = fb_filt[1] * fb_filt_ip[1]
    else
      fb_filt_op = vec(fb_filt)'*vec(fb_filt_ip);
    end
  end
  quantizer_op[1:end], SD[1:end]
end



  

