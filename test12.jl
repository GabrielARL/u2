using SignalAnalysis: length
using SignalAnalysis

tx=convert.(Int64, mseq(12))

for i = 1:length(tx)
    if(tx[i] == -1)
        tx[i] = 0
    end
end

testdata = zeros(Int64,1, 1024)
testdata1 = Main.FEC.encode(Code(63, 16, 11, 0o6331141367235453), testdata)












