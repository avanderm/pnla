module PNLA

    export
    #Types polyutil.jl
        PolySys,
        SymSys,
        MonomialOrder

    export
    #Types errors.jl
        OverDeterminedError,
        UnderDeterminedError

    export
    #Functions polyutil.jl
        plex,
        grlex,
        tdeg,
        reorder,
        sym2poly,
        poly2sym,
        poly2horner,
        evalsys

    export
    #Functions basics.jl
        degeqs,
        degmax
        
    export
    #Functions bounds.jl
        bezout

    include(joinpath("util", "polyutil.jl"))
    include("basics.jl")
    include("bounds.jl")
    include("errors.jl")

end #module
