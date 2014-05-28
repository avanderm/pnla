include("basics.jl")

function bezout(polysys::PolySys)
    # Number of variables must be equal to the number of equations
    @assert size(polysys.expn[1], 2) == length(polysys.expn)

    reduce(*, degeqs(polysys))
end
