module PolyUtil

using Base.Order

export
#Types
    PolySys,
    SymSys,
    MonomialOrder,
    NestedHorner

export
#Functions
    plex,
    grlex,
    tdeg,
    reorder,
    sym2poly,
    poly2sym,
    evalsys

const wp = Float64              # working precision
const ep = Int64                # exponent precision
const cp = Complex{Float64}     # calculation precision (variable values)

# Represents a polynomial system using coefficients and monomial exponents, variables are considered nameless
# from the numerical point of view
type PolySys
    coef::Array{Vector{wp},1}
    expn::Array{SparseMatrixCSC{ep,Int64},1}
end

# Represents a symbolic system of expressions in a number of user defined variables. Constructs a function
# implementing the behavior for easy evaluation
type SymSys
    expr::Expr
    vars::Expr

    f::Function

    function SymSys(expr::Expr, vars::Expr)
        f = eval(:($vars->$expr))
        new(expr,vars,f)
    end

    # Univariate case (Symbol is not a Expr)
    SymSys(expr::Expr, var::Symbol) = SymSys(expr, symbols(var))

    # Overloaded functions for text arguments
    SymSys(expr::Expr, vars::ASCIIString)        = SymSys(expr, parse(vars))
    SymSys(expr::ASCIIString, vars::ASCIIString) = SymSys(parse(expr), parse(vars))
end

# Represents a monomial order [Cox,2005]. A weight matrix m transforms the monomial order into a lexicographic
# ordering problem
type MonomialOrder
    m::SparseMatrixCSC{wp,Int64}
end

# Constructs the lexicographic order
function plex(ordering::Expr)
    return (vars::Expr)->plex_impl(vars, ordering)
end

plex(ordering::Symbol)      = plex(symbols(ordering))
plex(ordering::ASCIIString) = plex(parse(ordering))

# Implements the lexicographic order
function plex_impl(vars::Expr, ordering::Expr)
    @assert length(vars.args) == length(ordering.args)

    n::Int64 = length(ordering.args)
    if vars == ordering.args
        return MonomialOrder(speye(wp,n))
    else
        dict::Dict{Symbol,Int64} = Dict(vars.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering.args[i]] for i = 1:n ]

        return MonomialOrder(sparse(1:n, perm, ones(wp,n)))
    end
end

# Constructs the graded lexicographic order
function grlex(ordering::Expr)
    return (vars::Expr)->grlex_impl(vars, ordering)
end

grlex(ordering::Symbol)      = grlex(:($ordering,))
grlex(ordering::ASCIIString) = grlex(parse(ordering))

# Implements the graded lexicographic order
function grlex_impl(vars::Expr, ordering::Expr)
    @assert length(vars.args) == length(ordering.args)

    n::Int64 = length(ordering.args)
    if vars == ordering.args
        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,1:n-1],ones(wp,2*n-1)))
    else
        dict::Dict{Symbol,Int64} = Dict(vars.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering.args[i]] for i = 1:n-1 ]

        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,perm],ones(wp,2*n-1)))
    end
end

# Constructs the graded reverse lexicographic order
function tdeg(ordering::Expr)
    return (vars::Expr)->tdeg_impl(vars, ordering)
end

tdeg(ordering::Symbol)      = tdeg(:($ordering,))
tdeg(ordering::ASCIIString) = tdeg(parse(ordering))

# Implements the graded reverse lexicographic order
function tdeg_impl(vars::Expr, ordering::Expr)
    @assert length(vars.args) == length(ordering.args)

    n::Int64 = length(ordering.args)
    if vars == ordering.args
        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,n:-1:2],[ones(wp,n),-ones(wp,n-1)]))
    else
        dict::Dict{Symbol,Int64} = Dict(vars.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering.args[i]] for i = 2:n ]

        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,perm[end:-1:1]],[ones(wp,n),-ones(wp,n-1)]))
    end
end

# Reorder variables for a SymSys object (function needs to be reestablished for correct variable linking)
function reorder(sys::SymSys, vars::Expr)
    SymSys(sys.expr, vars)
end

# Builds a PolySys object using a SymSys object. Alongside construction, the equations are expanded,
# simplified, merged for equal powers with a user defined monomial ordering.
function sym2poly(symsys::SymSys, m::MonomialOrder)
    exprtype::Symbol = symsys.expr.head

    if exprtype == :tuple
        # System of s > 1 equations
        s = length(symsys.expr.args)

        coef = Array(Vector{wp},s)
        expn = Array(SparseMatrixCSC{ep,Int64},s)

        dict::Dict{Symbol,Integer} =
            Dict(convert(Array{Symbol,1}, symsys.vars.args), 1:length(symsys.vars.args))

        for i = 1:length(symsys.expr.args)
            (coef[i], expn[i]) = process_term(symsys.expr.args[i], sym2poly_dict(), dict)
            (expn[i], p) = monomial_sort(expn[i], m)
            (coef[i], expn[i]) = poly_reduce(coef[i][p], expn[i]) 
        end

        return PolySys(coef,expn)
    elseif exprtype == :call
        # System of s = 1 equations (embed in tuple and resend)
        return sym2poly(SymSys(Expr(:tuple,symsys.expr),symsys.vars), m)
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
        exp = spzeros(ep,1,length(vardict))
        exp[vardict[symbol]] = convert(ep,1)
        return [convert(wp,1.0)], exp 
    catch
        error("Did not recognize variable $(symbol), not a user defined symbol")
    end
end

# SYM2POLY
function sym2poly_number(number::Number, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    # Convert to working precision
    return [convert(wp, number)], spzeros(ep,1,length(vardict))
end

# SYM2POLY
function sym2poly_plus(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}) = process_term(terms[1], opdict, vardict)

    # Dynamic array growth
    for i = 2:length(terms)
        (coeft::Vector{wp}, expnt::SparseMatrixCSC{ep,Int64}) = process_term(terms[i], opdict, vardict)

        coef = vcat(coef, coeft)
        expn = vcat(expn, expnt)
    end

    return coef, expn
end

# SYM2POLY
function sym2poly_minus(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}) = process_term(terms[1], opdict, vardict)

    if length(terms) == 1
        # Expression of form -y
        return -coef, expn
    else
        # Expression of form x - y
        (coefr::Vector{wp}, expnr::SparseMatrixCSC{ep,Int64}) = process_term(terms[2], opdict, vardict)

        coef = vcat(coef, -coefr)
        expn = vcat(expn, expnr)
        return coef, expn
    end
end

# SYM2POLY
function sym2poly_product(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}) = process_term(terms[1], opdict, vardict)

    for i = 2:length(terms)
        (coeft::Vector{wp}, expnt::SparseMatrixCSC{ep,Int64}) = process_term(terms[i], opdict, vardict)

        # Coefficients can be determined by generalized outer product
        coef = kron(coeft, coef)
        # Exponents by matrix repetition and blockwise adding exponents of all monomials in the second term
        expn = vcat([ broadcast(+,expn,expnt[i,:]) for i in 1:size(expnt,1) ]...)
    end

    return coef, expn
end

# SYM2POLY
function sym2poly_power(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}) = process_term(terms[1], opdict, vardict)

    # Same principle as sym2poly_product but more efficiently implemented
    (coeft::Vector{wp}, expnt::SparseMatrixCSC{ep,Int64}) = (coef, expn)
    for i = 2:terms[2]
        coef = kron(coeft, coef)
        expn = vcat([ broadcast(+,expn,expnt[i,:]) for i in 1:size(expnt,1) ]...)
    end

    return coef, expn
end

# SYM2POLY
function sym2poly_divide(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}) = process_term(terms[1], opdict, vardict)

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

# Sorts monomials by transforming the order problem into a lexicographic setting using the monomial order
# weighting matrix [Cox,2005]
function monomial_sort(expn::SparseMatrixCSC{ep,Int64}, morder::MonomialOrder)
    # Transpose for SparseMatrixCSC is simple swap of two index vectors
    A = morder.m*expn'

    # Code for sortcols has been adopted from base library to obtain the permutation p
    r = 1:size(A,1)
    cols = [ sub(A,r,i) for i = 1:size(A,2) ]
    p = flipud(sortperm(cols, order=Lexicographic))

    return expn[p,:], p
end

# Merges coefficients of duplicate monomials in polynomial, assumes ordered by monomial_sort
function poly_reduce(coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64})
    n = length(coef)
    coefr::SparseMatrixCSC{wp,Int64} = spzeros(wp,n,1)

    current = 1
    coefr[current] = coef[current]

    for i = 2:n
        if expn[i,:] == expn[current,:]
            coefr[current] += coef[i]
        else
            current = i
            coefr[current] = coef[current]
        end
    end

    # Reduction of polynomial
    (Ti,Tj,Tz) = findnz(coefr)

    return Tz, expn[Ti,:]
end

# Builds a SymSys object using a PolySys object. Minus operators are supressed in sums by incorporating them
# into the coefficients of a monomial.
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
function poly2sym_eq(coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}, vars::Expr)
    nmons = length(coef)

    # Cannot assume eq to be of type Expression
    if nmons == 1
        eq = poly2sym_term(coef[1], expn[1,:], vars)
    else
        eq = Expr(:call, :+, [ poly2sym_term(coef[i], expn[i,:], vars)
            for i in 1:nmons ]...)
    end

    return eq
end

# A variant to poly2sym, returning a nested Horner system, better suited for evaluation
type NestedHorner
    function NestedHorner(polysys::PolySys, vars::Expr, m::MonomialOrder)
        @assert size(polysys.expn[1],2) == length(vars.args)

        s = size(polysys.expn,1)
        # Generate expressions
        expr::Expr
        if s > 1
            expr = Expr(:tuple,[ poly2sym_eq_nested(polysys.coef[i], polysys.expn[i], vars, m)
                for i = 1:s ]...)
        else
            expr = poly2sym_eq_nested(polysys.coef[1], polysys.expn[1], vars, m)
        end

        SymSys(expr, vars)
    end

    NestedHorner(polysys::PolySys, vars::Expr, morder::Function) = NestedHorner(polysys, vars, morder(vars))

    NestedHorner(polysys::PolySys, vars::Expr)  = NestedHorner(polysys, vars, tdeg(vars))
    NestedHorner(polysys::PolySys, var::Symbol) = NestedHorner(polysys, symbols(var), tdeg(var))

    NestedHorner(polysys::PolySys, vars::ASCIIString, m::MonomialOrder) =
        NestedHorner(polysys, parse(vars), m)
    NestedHorner(polysys::PolySys, vars::ASCIIString, morder::Function) =
        NestedHorner(polysys, parse(vars), morder(parse(vars)))
    NestedHorner(polysys::PolySys, vars::ASCIIString) = 
        NestedHorner(polysys, parse(vars))

    function NestedHorner(polysys::PolySys)
        nvars = size(polysys.expn[1],2)
    
        # Generate anonymous variables
        vars::Expr = Expr(:tuple,[ symbol("x$i") for i = 1:nvars ]...)

        NestedHorner(polysys, vars)
    end
end

# Processes a single equation using a variant of the nested Horner scheme [Carnicer, 1990]. We assume a graded
# monomial order.
function poly2sym_eq_nested(coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}, vars::Expr, m::MonomialOrder)
    r = length(coef)

    # Initial sort (note: graded property is more important than secondary sorting)
    (expns::SparseMatrixCSC{ep,Int64}, p::Array{ep,1}) = monomial_sort(expn, m)

    # Initial Horner nodes: leaf nodes storing a Number as expression
    hnodes::Array{HornerNode,1} = [ HornerNode(expns[i,:], coef[p[i]]) for i = 1:r ]

    # In each iteration we seek the smallest distance for each monomial of highest degree compared to all other
    # monomials and choose the best (pruning the tree stepwise top down)
    for i = 1:r-1
        dmax = sum(expns[1,:])

        hdeg = 1
        # Tracker for minimum value, offset index of minimal comparison exponent and index of the highest degree
        # monomial for which it occurred
        minc = (typemax(ep), 0, hdeg)
        while sum(expns[hdeg,:]) == dmax && hdeg != r-(i-1)
            # Compute distances for expns[hdeg,:] relative to all expns[hdeg+1:end,:]
            d::Array{ep,2} = sum(abs(expns[hdeg+1:r-(i-1),:] - repmat(expns[hdeg,:],r-i-hdeg+1)), 2)

            # Minimum, update current minimum if smaller
            mind = findmin(d)
            if mind[1] < minc[1]
                minc = tuple(mind..., hdeg)
            end

            hdeg += 1
        end

        # Replace Horner node hdeg by the merged node
        hnodes[minc[3]] = HornerNode(hnodes[minc[3]], hnodes[minc[2] + minc[3]], vars)

        # Update arrays
        expns[minc[3],:] = hnodes[minc[3]].expn
        expns = expns[[1:(minc[2] + minc[3] - 1); (minc[2] + minc[3] + 1):r-(i-1)],:]
        hnodes = hnodes[ [1:(minc[2] + minc[3] - 1); (minc[2] + minc[3] + 1):r-(i-1)] ]

        # Resort arrays for next iteration
        (expns, p) = monomial_sort(expns, m)
        hnodes = hnodes[p]
    end

    # One node left, need to reduce its exponent to [0 0 ... 0]
    try
        return Expr(:call, :*, hnodes[1].expr, poly2sym_monomial(hnodes[1].expn, vars))
    catch
        # Zero degree monomial occurence
        return hnodes[1].expr
    end
end

# Horner node, equal to the monomial coefficient if a leaf, merges two nodes if an intermediate node
type HornerNode
    expn::SparseMatrixCSC{ep,Int64}     # Exponent that the node in the Horner scheme represents
    expr::Any                           # Node value based upon predecessor nodes (Expr or Number)

    # Leaf node
    HornerNode(expn::SparseMatrixCSC{ep,Int64}, a::wp) = new(expn, a)

    # Merged node
    function HornerNode(left::HornerNode, right::HornerNode, vars::Expr)
        # Compute GCD of nodes
        gcd = sparse(min(left.expn, right.expn))

        exprl::Any = 0      # Number or Expression
        exprr::Any = 0

        try
            exprl = Expr(:call, :*, left.expr, poly2sym_monomial(sparse(left.expn - gcd), vars))
        catch
            exprl = left.expr
        end

        try
            exprr = Expr(:call, :*, right.expr, poly2sym_monomial(sparse(right.expn - gcd), vars))
        catch
            exprr = right.expr
        end

        new(gcd, Expr(:call, :+, exprl, exprr))
    end
end

# Processes a monomial term
function poly2sym_term(coef::wp, expn::SparseMatrixCSC{ep,Int64}, vars::Expr)
    try
        return Expr(:call, :*, coef, poly2sym_monomial(expn, vars))
    catch
        # Degree zero monomial
        return coef
    end
end

# Processes a leading monomial
function poly2sym_monomial(expn::SparseMatrixCSC{ep,Int64}, vars::Expr)
    # Helper function for powers, variables of power 1 added without explicit power
    pow = (sym::Symbol, n::Integer)->(
        if n == 1
            return sym
        else
            return Expr(:call, :^, sym, n)
        end)

    I,J,Z = findnz(expn)

    if length(I) == 0
        error("Zero monomial not represented as an expression")
    end

    if length(I) == 1
        # Univariate monomial
        return pow(vars.args[J[1]], Z[1])
    else
        # Variables of power 0 are excluded
        term::Expr = Expr(:call, :*, [ pow(vars.args[J[i]], Z[i]) for i in 1:length(I) ]...) 
    end

    return term
end

# Evaluates equation system in points stored as rows in a matrix
function evalsys{T<:Number}(sys::SymSys, val::Array{T,2})
    vcat([ collect(sys.f(convert(Array{cp,2}, val[i,:])...))' for i in 1:size(val,1) ]...)
end

# Evaluates equation system using array of vectors
function evalsys{T<:Number}(sys::SymSys, val::Array{Vector{T},1})
    vcat([ collect(sys.f(convert(Vector{cp},val[i])...)') for i in 1:length(val) ]...)
end

evalsys(sys::SymSys, val::Array{Any,1}) = evalsys(sys, convert(Array{Vector{cp},1}, val))

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

# Helper functions
function symbols(sym::Symbol)
    return :($sym,)
end

end #module
