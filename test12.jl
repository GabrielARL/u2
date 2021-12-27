using SignalAnalysis

include("FEC.jl")
include("BCH.jl")

import Main.FEC
import Main.BCH

tx=convert.(Int64, mseq(12))

for i = 1:length(tx)
    if(tx[i] == -1)
        tx[i] = 0
    end
end

data = BitVector(rand(Bool, 1040))
EncData=FEC.encode(BCH.Code(63, 16, 11, 0o6331141367235453), data)
sc = xor.(EncData, tx)
EncDataR = xor.(tx, sc)
dataR,_=FEC.decode(BCH.Code(63, 16, 11, 0o6331141367235453), EncDataR)

println(sum(data .-dataR))












