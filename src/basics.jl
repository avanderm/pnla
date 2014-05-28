include("global.jl")
include(joinpath("util", "polyutil.jl"))

# Computes the degree of each system equation
function degeqs(polysys::PolySys)
    da = @parallel [ maximum(sum(polysys.expn[i], 2)) for i=1:length(polysys.expn) ]
    convert(Array{ep,1}, da)
end

# Computes the maximum degree of all system equations
function degmax(polysys::PolySys)
    reduce(maximum, degeqs(polysys))
end
