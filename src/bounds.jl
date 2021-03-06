include("basics.jl")
include("errors.jl")

# Computes the number of projective solutions
function bezout(polysys::PolySys)
    # Number of variables must be equal to the number of equations
    s::Integer = length(polysys.expn)
    n::Integer = size(polysys.expn[1], 2)

    s > n? throw(OverDeterminedError(polysys)): s < n? throw(UnderDeterminedError(polysys)): reduce(*, degeqs(polysys))
end
