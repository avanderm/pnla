import Base.Order.Lexicographic

include(joinpath("..", "global.jl"))

typealias Var Union(Symbol, Expr)

# Represents a polynomial system using coefficients and monomial exponents, variables are considered nameless
# from the numerical point of view.
immutable type PolySys
    coef::Array{Vector{wp},1}
    expn::Array{SparseMatrixCSC{ep,Int64},1}
end

# Represents a symbolic system of expressions in a number of user defined variables. Constructs a function
# implementing the behavior for easy evaluation.
immutable type SymSys
    expr::Expr
    vars::Var
    func::Function

    function SymSys(expr::Expr, vars::Var)
        new(expr, varexpr(vars), eval(:($vars->$expr)))
    end

    # Overloaded functions for text arguments
    SymSys(expr::Expr,        vars::ASCIIString) = SymSys(expr, parse(vars))
    SymSys(expr::ASCIIString, vars::ASCIIString) = SymSys(parse(expr), vars)
end

# Represents a monomial order [Cox,2005]. A weight matrix m transforms the monomial order into a lexicographic
# ordering problem.
immutable type MonomialOrder
    m::SparseMatrixCSC{wp,Int64}
end

# Constructs the lexicographic order
function plex(ordering::Var)
    return (vars::Var)->plex_impl(vars, ordering)
end

plex(ordering::ASCIIString) = plex(parse(ordering))

# Implements the lexicographic order
function plex_impl(vars::Var, ordering::Var)
    vars_expr::Expr = varexpr(vars)
    ordering_expr::Expr = varexpr(ordering)
    @assert length(vars_expr.args) == length(ordering_expr.args)

    n::Int64 = length(ordering_expr.args)
    if vars == ordering_expr.args
        return MonomialOrder(speye(wp,n))
    else
        dict::Dict{Symbol,Int64} = Dict(vars_expr.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering_expr.args[i]] for i = 1:n ]

        return MonomialOrder(sparse(1:n, perm, ones(wp,n)))
    end
end

# Constructs the graded lexicographic order
function grlex(ordering::Var)
    return (vars::Var)->grlex_impl(vars, ordering)
end

grlex(ordering::ASCIIString) = grlex(parse(ordering))

# Implements the graded lexicographic order
function grlex_impl(vars::Var, ordering::Var)
    vars_expr::Expr = varexpr(vars)
    ordering_expr::Expr = varexpr(ordering)
    @assert length(vars_expr.args) == length(ordering_expr.args)

    n::Int64 = length(ordering_expr.args)
    if vars == ordering_expr.args
        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,1:n-1],ones(wp,2*n-1)))
    else
        dict::Dict{Symbol,Int64} = Dict(vars_expr.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering_expr.args[i]] for i = 1:n-1 ]

        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,perm],ones(wp,2*n-1)))
    end
end

# Constructs the graded reverse lexicographic order
function tdeg(ordering::Var)
    return (vars::Var)->tdeg_impl(vars, ordering)
end

tdeg(ordering::ASCIIString) = tdeg(parse(ordering))

# Implements the graded reverse lexicographic order
function tdeg_impl(vars::Var, ordering::Var)
    vars_expr::Expr = varexpr(vars)
    ordering_expr::Expr = varexpr(ordering)
    @assert length(vars_expr.args) == length(ordering_expr.args)

    n::Int64 = length(ordering_expr.args)
    if vars == ordering_expr.args
        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,n:-1:2],[ones(wp,n),-ones(wp,n-1)]))
    else
        dict::Dict{Symbol,Int64} = Dict(vars_expr.args, 1:n)
        perm::Vector{Int64} = [ dict[ordering_expr.args[i]] for i = 2:n ]

        return MonomialOrder(sparse([ones(Int64,n),2:n],[1:n,perm[end:-1:1]],[ones(wp,n),-ones(wp,n-1)]))
    end
end

# Reorder variables for a SymSys object (function needs to be reestablished for correct variable linking)
function reorder(sys::SymSys, vars::Var)
    SymSys(sys.expr, vars)
end

reorder(sys::SymSys, vars::ASCIIString) = reorder(sys, parse(vars))

# Builds a PolySys object using a SymSys object. Alongside construction, the equations are expanded,
# simplified, merged for equal powers with a user defined monomial ordering.
function sym2poly(symsys::SymSys, m::MonomialOrder)
    exprtype::Symbol = symsys.expr.head

    if exprtype == :tuple
        # System of s > 1 equations
        s::Integer = length(symsys.expr.args)

        coef::Array{Vector{wp},1} = Array(Vector{wp},s)
        expn::Array{SparseMatrixCSC{ep,Int64},1} = Array(SparseMatrixCSC{ep,Int64},s)

        dict::Dict{Symbol,Integer} =
            Dict(convert(Array{Symbol,1}, symsys.vars.args), 1:length(symsys.vars.args))

        for i = 1:length(symsys.expr.args)
            (coef[i], expn[i]) = process_term(symsys.expr.args[i], sym2poly_dict(), dict)
            (expn[i], p) = monomial_sort(expn[i], m)
            (coef[i], expn[i]) = poly_reduce(coef[i][p], expn[i]) 
        end

        return PolySys(coef, expn)
    elseif exprtype == :call
        # System of s = 1 equations (embed in tuple and resend)
        return sym2poly(SymSys(Expr(:tuple,symsys.expr),symsys.vars), m)
    else
        error("Invalid system type: $(symsys.expr.head)")
    end
end

# Monomial ordering occurs relative to the order stored by the SymSys object.
sym2poly(symsys::SymSys, mfun::Function = tdeg(symsys.vars)) = sym2poly(symsys, mfun(symsys.vars))

#SYM2POLY
function sym2poly_symbol(symbol::Symbol, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    # Coefficient 1 and unity vector
    try
        exp::SparseMatrixCSC{ep,Int64} = spzeros(ep,1,length(vardict))
        exp[vardict[symbol]] = convert(ep,1)
        return [ convert(wp,1.0) ], exp 
    catch
        error("Did not recognize variable $(symbol), not a user defined symbol")
    end
end

#SYM2POLY
function sym2poly_number(number::Number, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    # Convert to working precision
    return [ convert(wp, number) ], spzeros(ep,1,length(vardict))
end

#SYM2POLY
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

#SYM2POLY
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

#SYM2POLY
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

#SYM2POLY
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

#SYM2POLY
function sym2poly_divide(terms::Array{Any,1}, opdict::Dict{Symbol,Function}, vardict::Dict{Symbol,Integer})
    (coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}) = process_term(terms[1], opdict, vardict)

    # Divide all coefficients by what is assumed to be a Number
    try
        return (coef/eval(terms[2])), expn
    catch e
        error("Non valid operation in polynomial construction: division by $(e.var)")
    end
end

#SYM2POLY (circumvent Julia's not defined error)
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
function monomial_sort(expn::SparseMatrixCSC{ep,Int64}, mfun::MonomialOrder)
    # Transpose for SparseMatrixCSC is simple swap of two index vectors
    A::SparseMatrixCSC{ep,Int64} = mfun.m*expn'

    # Code for sortcols has been adopted from base library to obtain the permutation p
    r::UnitRange{Int64} = 1:size(A,1)
    cols = [ sub(A,r,i) for i = 1:size(A,2) ]
    p::Vector{ep} = flipud(sortperm(cols, order=Lexicographic))

    return expn[p,:], p
end

# Merges coefficients of duplicate monomials in polynomial, assumes ordered by monomial_sort
function poly_reduce(coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64})
    n::Integer = length(coef)
    coefr::SparseMatrixCSC{wp,Int64} = spzeros(wp,n,1)

    current::Integer = 1
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
    (Ti,Tj,Tz::Vector{wp}) = findnz(coefr)

    return Tz, expn[Ti,:]
end

# Builds a SymSys object using a PolySys object. Minus operators are supressed in sums by incorporating them
# into the coefficients of a monomial.
function poly2sym(polysys::PolySys, vars::Var, m::MonomialOrder)
    vars_expr::Expr = varexpr(vars)
    @assert size(polysys.expn[1],2) == length(vars_expr.args)

    s = size(polysys.expn,1)
    # Generate expressions
    expr::Expr
    if s > 1
        expr = Expr(:tuple,[ poly2sym_eq(polysys.coef[i], polysys.expn[i], vars_expr, m)
            for i = 1:s ]...)
    else
        expr = poly2sym_eq(polysys.coef[1], polysys.expn[1], vars_expr, m)
    end

    SymSys(expr,vars_expr)
end

# Function and symbol overloading
poly2sym(polysys::PolySys, vars::Var, mfun::Function = tdeg(vars)) = poly2sym(polysys, vars, mfun(vars))

# Text input overloading
poly2sym(polysys::PolySys, vars::ASCIIString, m::MonomialOrder) = poly2sym(polysys, parse(vars), m)
poly2sym(polysys::PolySys, vars::ASCIIString, mfun::Function = tdeg(parse(vars))) = poly2sym(polysys, vars, mfun(parse(vars)))

# Creates anonymous variables
function poly2sym(polysys::PolySys)
    n = size(polysys.expn[1],2)
    
    # Generate anonymous variables
    vars::Expr = Expr(:tuple,[ symbol("x$i") for i = 1:n ]...)

    poly2sym(polysys, vars)
end

# Processes a single equation, symbol or constant
function poly2sym_eq(coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}, vars::Expr, m::MonomialOrder)
    r = length(coef)

    # Initial sort
    (expns::SparseMatrixCSC{ep,Int64}, p::Array{ep,1}) = monomial_sort(expn, m)

    # Cannot assume eq to be of type Expression
    if r == 1
        eq = poly2sym_term(coef[p[1]], expns[1,:], vars)
    else
        eq = Expr(:call, :+, [ poly2sym_term(coef[p[i]], expns[i,:], vars)
            for i in 1:r ]...)
    end

    return eq
end

# A variant to poly2sym, returning a nested Horner system, better suited for evaluation
function poly2horner(polysys::PolySys, vars::Var, m::MonomialOrder)
    vars_expr::Expr = varexpr(vars)
    @assert size(polysys.expn[1],2) == length(vars_expr.args)

    s::Integer = size(polysys.expn,1)

    # Generate expressions
    expr::Expr
    if s > 1
        expr = Expr(:tuple,[ poly2sym_eq_nested(polysys.coef[i], polysys.expn[i], vars_expr, m)
            for i = 1:s ]...)
    else
        expr = poly2sym_eq_nested(polysys.coef[1], polysys.expn[1], vars_expr, m)
    end

    SymSys(expr, vars)
end

poly2horner(polysys::PolySys, vars::Var, mfun::Function = tdeg(vars)) = poly2horner(polysys, vars, mfun(vars))

poly2horner(polysys::PolySys, vars::ASCIIString, m::MonomialOrder) = poly2horner(polysys, parse(vars), m)
poly2horner(polysys::PolySys, vars::ASCIIString, mfun::Function = tdeg(parse(vars))) = poly2horner(polysys, vars, mfun(parse(vars)))

# Creates anonymous variables
function poly2horner(polysys::PolySys)
    n::Integer = size(polysys.expn[1],2)
    
    # Generate anonymous variables
    vars::Expr = Expr(:tuple,[ symbol("x$i") for i = 1:n ]...)

    poly2horner(polysys, vars)
end

# Processes a single equation using a variant of the nested Horner scheme [Carnicer, 1990]. We assume a graded
# monomial order.
function poly2sym_eq_nested(coef::Vector{wp}, expn::SparseMatrixCSC{ep,Int64}, vars::Expr, m::MonomialOrder)
    r::Integer = length(coef)

    # Initial sort (note: graded property is more important than secondary sorting)
    (expns::SparseMatrixCSC{ep,Int64}, p::Array{ep,1}) = monomial_sort(expn, m)

    # Initial Horner nodes: leaf nodes storing a Number as expression
    hnodes::Array{HornerNode,1} = [ HornerNode(expns[i,:], coef[p[i]]) for i = 1:r ]

    # In each iteration we seek the smallest distance for each monomial of highest degree compared to all other
    # monomials and choose the best (pruning the tree stepwise top down)
    for i = 1:r-1
        dmax::Integer = sum(expns[1,:])

        hdeg::Integer = 1
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
    expr::Union(Number, Expr)           # Node value based upon predecessor nodes (Expr or Number)

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

# Processes a monomial
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

# Evaluates equation system in points stored as columns in a matrix
function evalsys{T<:Number}(sys::SymSys, val::Array{T,2})
    @assert length(sys.vars.args) == size(val,1)
    valcp::Array{cp,2} = convert(Array{cp,2}, val)

    vcat([ collect(sys.func(valcp[:,i]...))' for i in 1:size(valcp,2) ]...)
end

# Evaluates equation system using 1-dim array of vectors
function evalsys(sys::SymSys, val::Array{Any,1})
    evalsys(sys, convert(Array{cp,2}, cat(2, val...)))
end

# Evaluates equations system using a variable number of vectors
evalsys{T<:Number}(sys::SymSys, x::Vector{T}...) = evalsys(sys, cat(2, x...))

# Special case for univariate polynomial evaluation (+ var. arg. fix)
evalsys(sys::SymSys,) = error("No evaluation points given")
evalsys(sys::SymSys, x::Number...) = evalsys(sys, convert(Array{cp,2}, collect(x)'))

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
function varexpr(vars::Var)
    if typeof(vars) == Expr
        return vars
    else
        return :($vars,)
    end
end
