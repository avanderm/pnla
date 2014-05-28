include("global.jl")
include(joinpath("util", "polyutil.jl"))

function degeqs(polysys::PolySys)
    da = @parallel [ maximum(sum(polysys.expn[i], 2)) for i=1:length(polysys.expn) ]
    convert(Array{ep,1}, da)
end

function degmax(polysys::PolySys)
    reduce(maximum, degeqs(polysys))
end
