module PolyUtil

using Base.Order

export
#Types
    PolySys,
    SymSys,
    MonomialOrder

export
#Functions
    plex,
    grlex,
    tdeg,
    reorder,
    sym2poly,
    poly2sym,
    evalsys

# Represents a polynomial system using coefficients and monomial exponents,
# variables are considered nameless from the numerical point of view
type PolySys
    coef::Array{Vector{Float64},1}
    expn::Array{Array{Int64,2},1}
end

# Represents a symbolic system of expressions in a number of user defined
# variables. Constructs a function implementing the behavior for easy evaluation
type SymSys
    expr::Expr
    vars::Expr

    f::Function

    function SymSys(expr::Expr, vars::Expr)
        f = eval(:($vars->$expr))
        new(expr,vars,f)
    end

    # Overloaded functions for text arguments
    function SymSys(expr::Expr, vars::ASCIIString)
        SymSys(expr, parse(vars))
    end

    function SymSys(expr::ASCIIString, vars::ASCIIString)
        SymSys(parse(expr), parse(vars))
    end
end

# Represents a monomial order [Cox,2005]. A weight matrix m transforms the
# monomial order into a lexicographic ordering problem
type MonomialOrder{T<:Number}
    m::SparseMatrixCSC{T,Int64}
end

# Represents the lexicographic order
function plex(ordering::Expr)
    return (vars::Expr)->plex_impl(vars, ordering)
end

function plex(ordering::ASCIIString)
    return plex(parse(ordering))
end

# Implements the lexicographic order
function plex_impl(vars::Expr, ordering::Expr)
    @assert length(vars.args) == length(ordering.args)

    n::Int64 = length(ordering.args)
    if vars == ordering.args
        return MonomialOrder(speye(n))
    else
        dict::Dict{Symbol,Int64} = Dict(vars.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering.args[i]] for i = 1:n ]

        return MonomialOrder(sparse(1:n, perm, ones(n)))
    end
end

# Represents the graded lexicographic order
function grlex(ordering::Expr)
    return (vars::Expr)->grlex_impl(vars, ordering)
end

function grlex(ordering::ASCIIString)
    return grlex(parse(ordering))
end

# Implements the graded lexicographic order
function grlex_impl(vars::Expr, ordering::Expr)
    @assert length(vars.args) == length(ordering.args)

    n::Int64 = length(ordering.args)
    if vars == ordering.args
        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,1:n-1],ones(2*n-1)))
    else
        dict::Dict{Symbol,Int64} = Dict(vars.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering.args[i]] for i = 1:n-1 ]

        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,perm],ones(2*n-1)))
    end
end

# Represents the graded reverse lexicographic order
function tdeg(ordering::Expr)
    return (vars::Expr)->tdeg_impl(vars, ordering)
end

function tdeg(ordering::ASCIIString)
    return tdeg(parse(ordering))
end

# Implements the graded reverse lexicographic order
function tdeg_impl(vars::Expr, ordering::Expr)
    @assert length(vars.args) == length(ordering.args)

    n::Int64 = length(ordering.args)
    if vars == ordering.args
        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,n:-1:2],[ones(n),-ones(n-1)]))
    else
        dict::Dict{Symbol,Int64} = Dict(vars.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering.args[i]] for i = 2:n ]

        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,perm[end:-1:1]],[ones(n),-ones(n-1)]))
    end
end


# Reorder variables for a SymSys object (function needs to be reestablished
# for correct variable linking)
function reorder(sys::SymSys, vars::Expr)
    SymSys(sys.expr, vars)
end

# Builds a PolySys object using a SymSys object. Alongside construction, the
# equations are expanded, simplified, merged for equal powers with a user
# defined monomial ordering.
function sym2poly(symsys::SymSys, m::MonomialOrder)
    exprtype::Symbol = symsys.expr.head

    if exprtype == :tuple
        # System of s > 1 equations
        s = length(symsys.expr.args)

        coef = Array(Vector{Float64},s)
        expn = Array(Array{Int64,2},s)

        dict::Dict{Symbol,Integer} =
            Dict(convert(Array{Symbol,1}, symsys.vars.args), 1:length(symsys.vars.args))

        for i = 1:length(symsys.expr.args)
            (coef[i], expn[i]) = poly_sort(poly_reduce(
                process_term(symsys.expr.args[i],sym2poly_dict(),dict)...)...,m)
        end

        return PolySys(coef,expn)
    elseif exprtype == :call
        # System of s = 1 equations (embed in tuple and resend)
        return sym2poly(SymSys(Expr(:tuple,symsys.expr),symsys.vars))
    else
        error("Invalid system type: $(symsys.expr.head)")
    end
end

# Monomial ordering occurs relative to the order stored by the SymSys object.
function sym2poly(symsys::SymSys, morder::Function)
    sym2poly(symsys, morder(symsys.vars))
end

# SYM2POLY
function sym2poly_symbol(symbol::Symbol, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    # Coefficient 1 and unity vector
    try
        exp = zeros(Int64,1,length(vardict))
        exp[vardict[symbol]] = 1
        return [1.0], exp 
    catch
        error("Did not recognize variable $(symbol), not a user defined symbol")
    end
end

# SYM2POLY
function sym2poly_number(number::Number, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    # Convert to double precision
    return [convert(Float64, number)], zeros(Int64,1,length(vardict))
end

# SYN2POLY
function sym2poly_plus(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{Float64}, expn::Array{Int64,2}) = process_term(terms[1], opdict, vardict)

    # Dynamic array growth
    for i = 2:length(terms)
        (coeft::Vector{Float64}, expnt::Array{Int64,2}) = process_term(terms[i], opdict, vardict)

        coef = vcat(coef, coeft)
        expn = vcat(expn, expnt)
    end

    return coef, expn
end

# SYM2POLY
function sym2poly_minus(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{Float64}, expn::Array{Int64,2}) = process_term(terms[1], opdict, vardict)

    if length(terms) == 1
        # Expression of form -y
        return -coef, expn
    else
        # Expression of form x - y
        (coefr::Vector{Float64}, expnr::Array{Int64,2}) = process_term(terms[2], opdict, vardict)

        coef = vcat(coef, -coefr)
        expn = vcat(expn, expnr)
        return coef, expn
    end
end

# SYM2POLY
function sym2poly_product(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{Float64}, expn::Array{Int64,2}) = process_term(terms[1], opdict, vardict)

    for i = 2:length(terms)
        (coeft::Vector{Float64}, expnt::Array{Int64,2}) = process_term(terms[i], opdict, vardict)

        # Coefficients can be determined by generalized outer product
        coef = kron(coeft, coef)
        # Exponents by matrix repetition and blockwise adding exponents of all monomials in
        # the second term
        expn = vcat([ broadcast(+,expn,expnt[i,:]) for i in 1:size(expnt,1) ]...)
    end

    return coef, expn
end

# SYM2POLY
function sym2poly_power(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{Float64}, expn::Array{Int64,2}) = process_term(terms[1], opdict, vardict)

    # Same principle as sym2poly_product but more efficiently implemented
    (coeft::Vector{Float64}, expnt::Array{Int64,2}) = (coef, expn)
    for i = 2:terms[2]
        coef = kron(coeft, coef)
        expn = vcat([ broadcast(+,expn,expnt[i,:]) for i in 1:size(expnt,1) ]...)
    end

    return coef, expn
end

# SYM2POLY
function sym2poly_divide(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{Float64}, expn::Array{Int64,2}) = process_term(terms[1], opdict, vardict)

    # Divide all coefficients by what is assumed to be a Number
    return (coef/eval(terms[2])), expn
end

# SYM2POLY (circumvent Julia's not defined error)
function sym2poly_dict()
    dict = {
        :+ => sym2poly_plus,
        :- => sym2poly_minus,
        :* => sym2poly_product,
        :^ => sym2poly_power,
        :/ => sym2poly_divide,
        :x => sym2poly_symbol,
        :c => sym2poly_number
    }

    return convert(Dict{Symbol,Function}, dict)
end

# Merges coefficients of duplicate monomials in polynomial
function poly_reduce(coef::Vector{Float64}, expn::Array{Int64,2})
    n = length(coef)
    coefc::Vector{Float64} = zeros(Float64,n)

    mask::Vector{Bool} = convert(Vector{Bool}, ones(n))

    for i = 1:n
        if ~mask[i]
            # Skip duplicates
            continue
        else
            ind::Vector{Int64} = filter(j->expn[j,:] == expn[i,:], i:n)

            # Update coefficient
            coefc[i] = sum(coef[ind])
            # Delete duplicates from mask
            if length(ind) > 1
                mask[ind[2:end]] = false
            end
        end
    end

    # Reduction of polynomial
    coefc = coefc[mask]
    expnc::Array{Int64,2} = expn[mask,:]

    return coefc, expnc
end

# Sorts the monomials by transforming the order problem into a lexicographic setting
# using the monomial order weighting matrix [Cox,2005]
function poly_sort(coef::Vector{Float64}, expn::Array{Int64,2}, morder::MonomialOrder)
    # Code for sortcols has been adopted from base library to obtain the permutation p
    A = morder.m*expn'

    r = 1:size(A,1)
    cols = [ sub(A,r,i) for i = 1:size(A,2) ]
    p = flipud(sortperm(cols, order=Lexicographic))

    return coef[p], expn[p,:]
end

# Builds a SymSys object using a PolySys object. Minus operators are supressed in sums by
# incorporating them into the coefficients of a monomial.
function poly2sym(polysys::PolySys, vars::Expr)
    @assert size(polysys.expn[1],2) == length(vars.args)

    nexpr = size(polysys.expn,1)
    # Generate expressions
    expr::Expr = Expr(:tuple,[ poly2sym_eq(polysys.coef[i], polysys.expn[i], vars)
        for i = 1:nexpr ]...)

    SymSys(expr,vars)
end

# Overloaded versions of poly2sym
function poly2sym(polysys::PolySys, vars::ASCIIString)
    poly2sym(polysys, parse(vars))
end

function poly2sym(polysys::PolySys)
    nvars = size(polysys.expn[1],2)
    
    # Generate anonymous variables
    vars::Expr = Expr(:tuple,[ symbol("x$i") for i = 1:nvars ]...)

    poly2sym(polysys, vars)
end

# Processes a single equation, symbol or constant
function poly2sym_eq(coef::Vector{Float64}, expn::Array{Int64,2}, vars::Expr)
    nmons = length(coef)

    # Cannot assume eq to be of type Expression
    if nmons == 1
        eq = poly2sym_monomial(coef[1], expn[1,:], vars)
    else
        eq = Expr(:call, :+, [ poly2sym_monomial(coef[i], expn[i,:], vars)
            for i in 1:nmons ]...)
    end

    return eq
end

# Processes a monomial
function poly2sym_monomial(coef::Float64, expn::Array{Int64,2}, vars::Expr)
    if findnext(expn,1) == 0
        # Constant case
        return coef
    else
        # Variables of power 0 are excluded, variables of power 1 added without explicit power
        monomial::Expr = Expr(:call, :*, coef, [
            ((sym::Symbol, n::Integer)->(
            if n == 1
                return sym
            else
                return Expr(:call, :^, sym, n)
            end))(vars.args[i], expn[i]) for i in find(expn) ]...) 

        return monomial
    end
end

# Evaluates equation system in points stores as row in a matrix
function evalsys(sys::SymSys, val::Array{Float64,2})
    [ sys.f(val[i,:]...) for i in 1:size(val,1) ]
end

function evalsys(sys::SymSys, val::Array{Any,2})
    evalsys(sys, convert(Array{Float64,2}, val))
end

# Evaluates equation system using array of vectors (double precision)
function evalsys(sys::SymSys, val::Array{Vector{Float64},1})
    [ sys.f(val[i]...) for i in 1:length(val) ]
end

function evalsys(sys::SymSys, val::Array{Any,1})
    evalsys(sys, convert(Array{Vector{Float64},1}, val))
end

# Processes term of type Expression and determines first operator
function process_term(eq::Expr, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    # Test if expression (no nested tuples, assignment, ...)
    if eq.head == :call
        # First operator
        if haskey(opdict, eq.args[1])
            return opdict[eq.args[1]](eq.args[2:end], opdict, vardict)
        else
            # Not recognized operator, assume Number and attempt evaluation
            return process_term(eval(eq))
        end
    else
        error("Invalid equation: $(eq.head)")
    end
end

# Processes term of type Symbol
function process_term(symbol::Symbol, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    if haskey(opdict, :x)
        return opdict[:x](symbol, opdict, vardict)
    end
end

# Processes term of type Number
function process_term(number::Number, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    if haskey(opdict, :c)
        return opdict[:c](number, opdict, vardict)
    end
end

end #module
