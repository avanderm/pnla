import Base.showerror

include(joinpath("util", "polyutil.jl"))

type OverDeterminedError <: Exception
    polysys::PolySys
end

showerror(io::IO, e::OverDeterminedError) = print(io,
    "overdetermined system s > n: $(length(e.polysys.expn)) equations for $(size(e.polysys.expn[1], 2)) variables")

type UnderDeterminedError <: Exception
    polysys::PolySys
end

showerror(io::IO, e::UnderDeterminedError) = print(io,
    "underdetermined system s < n: $(length(e.polysys.expn)) equations for $(size(e.polysys.expn[1], 2)) variables")
