moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

using ToeplitzMatrices

function makevb(rbb,d) 
#y = signal(repeat(mseq(12); inner=12) .* cw(-1000.0, length(mseq(12))*12/6000, 6000.0), 6000.0)
#yb = y .* cw(1000.0, length(mseq(12))*12/6000, 6000.0)  

yb = signal(repeat(mseq(12); inner=12), 6000)

rbb = resample(rbb,1000/999)
rbb = rbb./norm(rbb,4)
c = MMSE_eqlz(rbb, samples(yb), 12*d)
h = (1/c) 
h = h/norm(h,2)
h = h'
vn = zeros(ComplexF64, length(rbb), 1)
for i =length(h):length(rbb)
vn[i] = rbb[i] - vec(h[1:end])'*reverse(rbb[i-length(h)+1:1:i])
end
vn
end

function adaptRLSTrainPLL(s, bb, N_ff, N_fb, λ, modulation)
    ff_filt_len = N_ff
    fb_filt_len = N_fb
    Δ = 0.01
    training_len = length(s)
    ff_filt = zeros(ComplexF64, ff_filt_len)
    fb_filt = zeros(ComplexF64, fb_filt_len)
    ff_fb_filt = zeros(ComplexF64, fb_filt_len + ff_filt_len)
    ff_filt_ip = zeros(ComplexF64, ff_filt_len)
    fb_filt_ip = zeros(ComplexF64, fb_filt_len)
    SD = zeros(ComplexF64, training_len-ff_filt_len+1)
    error = zeros(ComplexF64, training_len-ff_filt_len+1)
    ff_filt_op = zeros(ComplexF64, training_len-ff_filt_len+1)
    fb_filt_op = zeros(ComplexF64, training_len-ff_filt_len+1)
    quantizer_op = zeros(ComplexF64, training_len-ff_filt_len+1)
    DO = zeros(Float64, training_len-ff_filt_len+1)
    ϕ = zeros(Float64, training_len-ff_filt_len+1)
    println(training_len)
    P = Δ * I(ff_filt_len+fb_filt_len)
    P′ = Δ * I(ff_filt_len+fb_filt_len)
    for i1=fb_filt_len+1:training_len-(fb_filt_len+ff_filt_len)
        ff_filt = ff_fb_filt[1:ff_filt_len]
        fb_filt = ff_fb_filt[(ff_filt_len+1):end]
        #println(i1)
        ff_filt_ip = s[(i1+ff_filt_len-1):-1:i1] * cis(mod(ϕ[i1-1], 2π))
        ff_filt_op = ((ff_filt)') * ff_filt_ip
        #println( (i1, fb_filt_len, i1-fb_filt_len) )
        
        fb_filt_ip = bb[(i1-1):-1:(i1-fb_filt_len)]
        fb_filt_op = ((fb_filt)') * fb_filt_ip
        combined_ip = [ff_filt_ip; fb_filt_ip]
        ff_and_fb = ff_filt_op + fb_filt_op
        SD[i1] = ff_and_fb
        temp1 = convert(Float64, sign(real(ff_and_fb)));
        if(modulation == 2)
            temp2 = convert(Float64, sign(imag(ff_and_fb)));
        else
            temp2 = 0
        end
        quantizer_op[i1] = temp1 + 1im*(temp2);
        error[i1] = bb[i1] - ff_and_fb
        nume1 = P * combined_ip
        denom1 = λ + (combined_ip)'*nume1
        K = nume1/ denom1
        
        P′ = K*(combined_ip)'*P
        
        #println(size(P - P′))
        
        P = (P - P′)/λ

        DO[i1] = imag(log(ff_filt_op*conj(-error[i1]+ff_filt_op)))
        ϕ[i1] = ϕ[i1-1] - ϕ[i1-2] + 0.001*DO[i1] - 0.0001*DO[i1-1]

        ff_fb_filt = ff_fb_filt + K*(error[i1])'
    end
    println(norm(error[1000:end],2))
    error=moving_average(error, 1000)
    quantizer_op, error,ff_filt, fb_filt, SD;
  end
  function adaptRLSTrain2(s, bb, N_ff, N_fb, λ, modulation)
    ff_filt_len = N_ff
    fb_filt_len = N_fb
    Δ = 0.01
    training_len = min(length(s),length(bb))
    ff_filt = zeros(ComplexF64, ff_filt_len)
    fb_filt = zeros(ComplexF64, fb_filt_len)
    ff_fb_filt = zeros(ComplexF64, fb_filt_len + ff_filt_len)
    ff_filt_ip = zeros(ComplexF64, ff_filt_len)
    fb_filt_ip = zeros(ComplexF64, fb_filt_len)
    SD = zeros(ComplexF64, training_len-ff_filt_len+1)
    error = zeros(ComplexF64, training_len-ff_filt_len+1)
    ff_filt_op = zeros(ComplexF64, training_len-ff_filt_len+1)
    fb_filt_op = zeros(ComplexF64, training_len-ff_filt_len+1)
    quantizer_op = zeros(ComplexF64, training_len-ff_filt_len+1)
    #println(training_len)
    P = Δ * I(ff_filt_len+fb_filt_len)
    P′ = Δ * I(ff_filt_len+fb_filt_len)

    fbops = zeros(ComplexF64, training_len-(ff_filt_len)+1)
    ffops = zeros(ComplexF64, training_len-(ff_filt_len)+1)

    DO = zeros(Float64, training_len-ff_filt_len+1)
    ϕ = zeros(Float64, training_len-ff_filt_len+1)

    for i1=fb_filt_len+1:training_len-(fb_filt_len+ff_filt_len) 
        ff_filt = ff_fb_filt[1:ff_filt_len]
        fb_filt = ff_fb_filt[(ff_filt_len+1):end]
        #println(i1)
        ff_filt_ip = s[(i1+ff_filt_len-1):-1:i1] #* cis(-mod(ϕ[i1-1], 2π))
        ff_filt_op = ((ff_filt)') * ff_filt_ip
        #println( (i1, fb_filt_len, i1-fb_filt_len) )
        fb_filt_ip = bb[(i1-1):-1:(i1-fb_filt_len)] 
        fb_filt_op = ((fb_filt)') * fb_filt_ip
        combined_ip = [ff_filt_ip; fb_filt_ip]
        ff_and_fb = ff_filt_op + fb_filt_op
        error[i1] = bb[i1] - ff_and_fb
        SD[i1] = (ff_and_fb)
        fbops[i1] = fb_filt_op
        ffops[i1] = ff_filt_op
        temp1 = convert(Float64, sign(real(ff_and_fb)));
        if(modulation == 2)
            temp2 = convert(Float64, sign(imag(ff_and_fb)));
        else
            temp2 = 0
        end
        quantizer_op[i1] = temp1 + 1im*(temp2);
       
        nume1 = P * combined_ip
        denom1 = λ + (combined_ip)'*nume1
        K = nume1/ denom1
        
        P′ = K*(combined_ip)'*P
        DO[i1] = imag(ff_filt_op*(conj(quantizer_op[i1]+fb_filt_op)))
        ϕ[i1] = ϕ[i1-1] + 0.001*DO[i1] - 0.0001*DO[i1-1]
        #println(size(P - P′))
        
        P = (P - P′)/λ
        ff_fb_filt = ff_fb_filt + K*(error[i1])'
    end
    #println(norm(error[1000:end],2))
    error=moving_average(error, 1000)
    quantizer_op, error,ff_filt, fb_filt, SD, ffops, fbops, ff_fb_filt;
  end

function adaptRLSTrain(s, bb, N_ff, N_fb, λ, modulation)
  ff_filt_len = N_ff
  fb_filt_len = N_fb
  Δ = 0.01
  training_len = length(s)
  ff_filt = zeros(ComplexF64, ff_filt_len)
  fb_filt = zeros(ComplexF64, fb_filt_len)
  ff_fb_filt = zeros(ComplexF64, fb_filt_len + ff_filt_len)
  ff_filt_ip = zeros(ComplexF64, ff_filt_len)
  fb_filt_ip = zeros(ComplexF64, fb_filt_len)
  SD = zeros(ComplexF64, training_len-ff_filt_len+1)
  error = zeros(ComplexF64, training_len-ff_filt_len+1)
  ff_filt_op = zeros(ComplexF64, training_len-ff_filt_len+1)
  fb_filt_op = zeros(ComplexF64, training_len-ff_filt_len+1)
  quantizer_op = zeros(ComplexF64, training_len-ff_filt_len+1)
  println(training_len)
  P = Δ * I(ff_filt_len+fb_filt_len)
  P′ = Δ * I(ff_filt_len+fb_filt_len)
  for i1=fb_filt_len+1:training_len-(fb_filt_len+ff_filt_len)
      ff_filt = ff_fb_filt[1:ff_filt_len]
      fb_filt = ff_fb_filt[(ff_filt_len+1):end]
      #println(i1)
      ff_filt_ip = s[(i1+ff_filt_len-1):-1:i1]
      ff_filt_op = ((ff_filt)') * ff_filt_ip
      #println( (i1, fb_filt_len, i1-fb_filt_len) )
      
      fb_filt_ip = bb[(i1-1):-1:(i1-fb_filt_len)]
      fb_filt_op = ((fb_filt)') * fb_filt_ip
      combined_ip = [ff_filt_ip; fb_filt_ip]
      ff_and_fb = ff_filt_op + fb_filt_op
      SD[i1] = ff_and_fb
      temp1 = convert(Float64, sign(real(ff_and_fb)));
      if(modulation == 2)
          temp2 = convert(Float64, sign(imag(ff_and_fb)));
      else
          temp2 = 0
      end
      quantizer_op[i1] = temp1 + 1im*(temp2);
      error[i1] = bb[i1] - ff_and_fb
      nume1 = P * combined_ip
      denom1 = λ + (combined_ip)'*nume1
      K = nume1/ denom1
      
      P′ = K*(combined_ip)'*P
      
      #println(size(P - P′))
      
      P = (P - P′)/λ
      ff_fb_filt = ff_fb_filt + K*(error[i1])'
  end
  println(norm(error[1000:end],2))
  error=moving_average(error, 1000)
  quantizer_op, error,ff_filt, fb_filt, SD;
end

function adaptDD(s, ff_filt ,fb_filt, modulation)
  ff_filt_len = length(ff_filt)
  fb_filt_len = length(fb_filt)
  data_len = length(s)
  ff_filt_ip = zeros(ComplexF64, ff_filt_len)
  fb_filt_ip = zeros(ComplexF64, fb_filt_len)
  quantizer_op = zeros(ComplexF64, data_len-ff_filt_len+1)
  fb_filt_op = 0;
  quantizer_op = Array{ComplexF64}(undef,(length(s)-length(ff_filt)))
  SD = zeros(ComplexF64, data_len-ff_filt_len+1)
  fbops = zeros(ComplexF64, data_len-(ff_filt_len)+1)
  ffops = zeros(ComplexF64, data_len-(ff_filt_len)+1)
  DO = zeros(Float64, data_len-ff_filt_len+1)
  ϕ = zeros(Float64, data_len-ff_filt_len+1)

  for i1=1+1:data_len-(ff_filt_len)
      ff_filt_ip = s[(i1+ff_filt_len-1):-1:i1] #* 
      ff_filt_op = ((ff_filt)') * ff_filt_ip 
      ff_filt_op = ff_filt_op #* cis(-mod(ϕ[i1-1], 2π))
      ff_and_fb = ff_filt_op + fb_filt_op
      #ff_and_fb = ff_and_fb + abs.(imag(ff_and_fb))*im
      SD[i1] = ff_and_fb
      fbops[i1] = fb_filt_op
      ffops[i1] = ff_filt_op
      temp1 = convert(Float64, sign(real(ff_and_fb)));
      if(modulation == 2)
          temp2 = convert(Float64, sign(imag(ff_and_fb)));
          quantizer_op[i1] = temp1 + 1im*(temp2);
      else
          temp2 = 0
          quantizer_op[i1] = temp1
      end
      
      #error[i1] = round(ff_and_fb - training_seq[i1]; sigdigits = 32)
      fb_filt_ip[2:end]=fb_filt_ip[1:end-1];
      fb_filt_ip[1] = quantizer_op[i1];#ff_and_fb
      fb_filt_op = (fb_filt)' * fb_filt_ip


      DO[i1] = imag(ff_filt_op*(conj(quantizer_op[i1]+fb_filt_op )))
      ϕ[i1] = ϕ[i1-1] + 0.001*DO[i1] - 0.0001*DO[i1-1]


      #println((i1,fb_filt_op))
      #if(mod(i1,10)==0)
      #display(plot!(real.(fb_filt_ip)))
      #display(plot(imag.(fb_filt_ip))) 
      #end
  end
  quantizer_op, SD, fbops, ffops;
end

function chRlsTrain(s, bb, N_ff, λ, modulation)
  ff_filt_len = N_ff
  Δ = 0.001
  training_len = min(length(s),length(bb))
  ff_filt = zeros(ComplexF64, ff_filt_len)
  ff_filt_ip = zeros(ComplexF64, ff_filt_len)
  SD = zeros(ComplexF64, training_len-ff_filt_len+1)
  error = zeros(ComplexF64, training_len-ff_filt_len+1)
  ff_filt_op = zeros(ComplexF64, training_len-ff_filt_len+1)
  quantizer_op = zeros(ComplexF64, training_len-ff_filt_len+1)
  P = Δ * I(ff_filt_len)
  P′ = Δ * I(ff_filt_len)
  ffops = zeros(ComplexF64, training_len-(ff_filt_len)+1)

  ff = ComplexF64[]

  for i1=1:training_len-(+ff_filt_len) 
    ff_filt_ip = s[(i1+ff_filt_len-1):-1:i1] #* cis(-mod(ϕ[i1-1], 2π))
    ff_filt_op = ((ff_filt)') * ff_filt_ip
    error[i1] = bb[i1] - ff_filt_op
    SD[i1] = (ff_filt_op)
    ffops[i1] = ff_filt_op
    temp1 = convert(Float64, sign(real(ff_filt_op)));
    if(modulation == 2)
      temp2 = convert(Float64, sign(imag(ff_filt_op)));
    else
      temp2 = 0
    end
      quantizer_op[i1] = temp1 + 1im*(temp2)
      nume1 = P * ff_filt_ip
      denom1 = λ + (ff_filt_ip)'*nume1
      K = nume1/ denom1
      P′ = K*(ff_filt_ip)'*P
      P = (P - P′)/λ
      ff_filt = ff_filt + K*(error[i1])'

      for i=1:length(K)
        push!(ff,K[i])
      end
    end
    error=moving_average(error, 1000)
    quantizer_op, error,ff_filt, SD, ffops, ff;
end




function adaptLMSTrain(s, bb, N_ff, N_fb, δ, modulation)
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
      ff_filt_op = transpose(ff_filt)*(ff_filt_ip)
      ff_and_fb = ff_filt_op-fb_filt_op;
      err[i1] = (ff_and_fb - bb[i1])
      temp1 = convert(Float64, sign(real(ff_and_fb)));
      if(modulation == 2)
          temp2 = convert(Float64, sign(imag(ff_and_fb)));
          quantizer_op[i1] = temp1 + 1im*(temp2);
      else
          temp2 = 0
          quantizer_op[i1] = temp1
      end

      ff_filt=ff_filt-δ*err[i1]*conj.(ff_filt_ip);
      fb_filt=fb_filt+δ*err[i1]*conj.(fb_filt_ip);
      fb_filt_ip[2:end]=fb_filt_ip[1:end-1];
      fb_filt_ip[1] = quantizer_op[i1];#(ff_and_fb)#
      fb_filt_op = transpose(fb_filt)*(fb_filt_ip);
    end
    #println(norm(err[1000:end],2))
    #err=moving_average(err, 1000)
    quantizer_op, ff_filt, fb_filt, err;
  end

  function downsymbol(yb)
    temp = ComplexF64[]
    for i = 12:12:length(yb)
      push!(temp,convert.(ComplexF64,mean(yb[i-11:i])))
    end
    temp
  end

  function estimateD( s, bb )
    h=mfilter(bb,s)
    h
  end

  function MMSE_eqlz(s, MSEQ, M)
    ϕff = xcorr(samples(s), samples(s))
    ϕfd = xcorr(samples(MSEQ), samples(s))
    N = length(MSEQ)
    𝑅ff = ϕff[N:N+M-1]
    𝑅 = Toeplitz(𝑅ff, vec(𝑅ff))
    P = ϕfd[N:N+M-1]
    h = inv(collect(𝑅)) * (P)
    h
  end

  function DA(s, ff_filt , modulation)
    ff_filt_len = length(ff_filt)
    data_len = length(s)
    ff_filt_ip = zeros(ComplexF64, ff_filt_len)
    quantizer_op = zeros(ComplexF64, data_len-ff_filt_len+1)
    quantizer_op = Array{ComplexF64}(undef,(length(s)-length(ff_filt)))
    SD = zeros(ComplexF64, data_len-ff_filt_len+1)
  
    for i1=1+ff_filt_len:1:data_len-(ff_filt_len)
        ff_filt_ip =  s[(i1+ff_filt_len-1):-1:i1]  
        ff_filt_op = (ff_filt)' * ff_filt_ip 
        ff_and_fb = ff_filt_op 
        #ff_and_fb = ff_and_fb + abs.(imag(ff_and_fb))*im
        SD[i1] = ff_and_fb
        temp1 = convert(Float64, sign(real(ff_and_fb)));
        if(modulation == 2)
            temp2 = convert(Float64, sign(imag(ff_and_fb)));
            quantizer_op[i1] = temp1 + 1im*(temp2);
        else
            temp2 = 0
            quantizer_op[i1] = temp1
        end
    end
    quantizer_op, SD;
  end