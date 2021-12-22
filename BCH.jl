module BCH

using GaloisFields
using LinearAlgebra
using ..FEC

struct Code <: FEC.Code
  n::Int
  k::Int
  t::Int
  poly::BigInt
end

Base.show(io::IO, bch::Code) = print(io, "BCH($(bch.k)/$(bch.n))")

const BCH_63_16 = Code(63, 16, 11, 0o6331141367235453)
const BCH_31_6 = Code(31, 6, 7, 0o313365047)
const BCH_31_11 = Code(31, 11, 5, 0o5423325)

const list = [
  BCH_63_16,
  BCH_31_6,
  BCH_31_11
]

const GF2 = GaloisField(2)

function rref!(A::Matrix{GF2})
  n = min(size(A)...)
  i = 0
  k = 1
  while k ≤ n && i < size(A,1)
    i += 1
    if iszero(A[i,k])
      for j ∈ i+1:size(A,1)
        if isone(A[j,k])
          tmp = A[i,:]
          A[i,:] .= A[j,:]
          A[j,:] .= tmp
          break
        end
      end
      isone(A[i,k]) || continue
    end
    for j ∈ 1:size(A,1)
      if i != j && isone(A[j,k])
        A[j,:] .+= A[i,:]
      end
    end
    k += 1
  end
end

function generator(bch::Code)
  G = zeros(GF2, bch.k, bch.n)
  p = convert(Vector{GF2}, Vector{UInt8}(string(bch.poly; base=2)) .- 0x30)
  for i ∈ 1:bch.k
    G[bch.k-i+1,bch.n-length(p)-i+2:bch.n-i+1] .= p
  end
  rref!(G)
  G
end

function paritycheck(bch::Code)
  G = generator(bch)
  hcat(transpose(G[:,bch.k+1:end]), I(bch.n - bch.k) * GF2(1))
end

checkblock(bch::Code, bits) = mapreduce(isone, +, paritycheck(bch) * bits) == 0
encodeblock(bch::Code, bits) = map(isone, transpose(generator(bch)) * bits)

function decodeblock(bch::Code, r)
  r̄ = harddecision.(r)
  checkblock(bch, r̄) && return r̄
  G = generator(bch)
  λ₁ = sortperm(-abs.(r))
  y = r[λ₁]
  G′ = G[:,λ₁]
  G′′ = G′[:,1]
  λ₂ = [1]
  for i ∈ 2:size(G′,2)
    size(G′′,2) ≥ bch.k && break
    T = hcat(G′′, G′[:,i])
    rref!(T)
    if mapreduce(isone, +, T[:,end]) > 0
      G′′ = hcat(G′′, G′[:,i])
      push!(λ₂, i)
    end
  end
  for i ∈ 2:size(G′,2)
    if i ∉ λ₂
      G′′ = hcat(G′′, G′[:,i])
      push!(λ₂, i)
    end
  end
  z = y[λ₂]
  G₁ = deepcopy(G′′)
  rref!(G₁)
  a = map(>(0), z[1:bch.k])
  ā = map(isone, transpose(G₁) * a)         # order-0 estimate
  score(b) = -sum((z[i] > 0) != ā[i] ? abs(z[i]) : 0.0 for i ∈ 1:length(z))
  best = (ā, score(ā))
  for i ∈ 1:bch.k
    for j ∈ i:bch.k
      for k ∈ j:bch.k
        i != j && j == k && continue
        ã = deepcopy(a)
        ã[i] ⊻= 1
        j != i && (ã[j] ⊻= 1)
        k != j && (ã[k] ⊻= 1)
        ā = map(isone, transpose(G₁) * ã)   # order-3 estimate
        s̄ = score(ā)
        s̄ > best[2] && (best = (ā, s̄))
      end
    end
  end
  best[1][invperm(λ₂)][invperm(λ₁)]
end

### FEC interface functions

function FEC.create(::Type{Code}, k, n)
  best = missing
  wasted = n - k
  for bch ∈ list
    blks = cld(k, blocklength(bch))
    nbitsout = codewordlength(bch) * blks
    if nbitsout ≤ n
      if n - nbitsout < wasted
        wasted = n - nbitsout
        best = bch
      end
    end
  end
  best
end

FEC.coderate(bch::Code) = bch.k // bch.n
FEC.blocklength(bch::Code) = bch.k
FEC.codewordlength(bch::Code) = bch.n

function FEC.encode(bch::Code, bits; interleave=true)
  length(bits) % bch.k == 0 || throw(ArgumentError("Data length not a multiple of $(bch.k)"))
  b = [encodeblock(bch, b1) for b1 ∈ Iterators.partition(bits, bch.k)]
  interleave ? (hcat(b...) |> transpose |> vec) : vcat(b...)
end

function FEC.decode(bch::Code, metrics; interleave=true)
  length(metrics) % bch.n == 0 || throw(ArgumentError("Data length not a multiple of $(bch.n)"))
  inbits = harddecision.(metrics)
  interleave && (metrics = reshape(metrics, :, bch.n) |> transpose |> vec)
  bits = vcat([decodeblock(bch, m1)[1:bch.k] for m1 ∈ Iterators.partition(metrics, bch.n)]...)
  rebits = encode(bch, bits; interleave)
  errs = count(inbits .⊻ rebits)
  bits, errs
end

end # BCH
