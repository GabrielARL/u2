module FEC

export coderate, blocklength, codewordlength, encode, decode
export pad!, harddecision, bytes2bits, bits2bytes

### parent type for all FEC codes

"""
Forward Error Correction Code.
"""
abstract type Code end

### interface functions

"""
    create(::Type{<:FEC.Code}, k, n)

Create a code which takes in `k` bits and generates at maximum `n` bits.
The effective code rate for such a code is no less than `k/n`. Returns
`missing` if no code matching the given specification is available.
"""
function create end

"""
    coderate(fec::FEC.Code)

Get the code rate for a `fec` code. The code rate is usually returned as a
`Rational` number, but any `Real` number is acceptable.
"""
function coderate end

"""
    blocklength(fec::FEC.Code)

Get the length in bits for a uncoded block. For convolutional codes, block length
is returned as 1.
"""
function blocklength end

"""
    codewordlength(fec::FEC.Code)

Get the length in bits for a coded block (codeword). For convolutional codes, codeword
length is returned as 1.
"""
function codewordlength end

"""
    encode(fec::FEC.Code, bits)

Encode `bits` (typically given as a `BitVector`) using the `fec` code and return
a vector of encoded bits.
"""
function encode end

"""
    decode(fec::FEC.Code, metrics)

Decode `metrics` (typically given as a `AbstractVector{<:Real}`) using the `fec` code
and return a vector of decoded bits and the number of errors corrected. If decoding
fails, `(nothing, nothing)` may be returned.

Metrics are positive for estimated `1` bits and negative for estimated `0` bits.
Examples of good metrics include log-likelihood ratios.
"""
function decode end

### utility functions

"""
    pad!(bits, fec::FEC.Code)

Pad bits with `0` to be an integral multiple of block length for the `fec` code.
"""
function pad!(bits, fec::Code)
  k = blocklength(fec)
  length(bits) % k == 0 && return bits
  append!(bits, zeros(eltype(bits), k - (length(bits) % k)))
end

"""
    harddecision(metric)

Hard decision decoding of `metric` into a bit. A metric is positive for an estimated `1`
bit and negative for an estimated `0` bit. Examples of a good metric includes log-likelihood ratio.
"""
harddecision(metric) = metric > zero(eltype(metric))

"""
    byte2bits(b)

Convert byte `b` to a vector of bits.
"""
byte2bits(b) = BitVector((b & (0x80 >> n)) != 0 for n ∈ 0:7)

"""
    bits2bytes(bits)

Convert bit sequence in vector `bits` to a byte.
"""
bits2byte(b8) = UInt8(mapreduce(x -> x[2] << (8-x[1]), |, enumerate(b8)))

"""
    bytes2bits(b)

Convert byte sequence in vector `b` to a vector of bits.
"""
bytes2bits(bytes) = vcat([byte2bits(b) for b ∈ bytes]...)

"""
    bits2bytes(bits)

Convert bit sequence in vector `bits` to a vector of bytes.
"""
bits2bytes(bits) = [bits2byte(b8) for b8 ∈ Iterators.partition(bits, 8)]

end # FEC
