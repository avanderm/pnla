pnla
====

Polynomial Numerical Linear Algebra package for Julia

# Introduction

This package is in development. For a great PNLA package in a MATLAB/Oktave environment, check out Kim Batselier's code project [here](https://github.com/kbatseli/PNLA_MATLAB_OCTAVE).

## Constructing systems

In order to start working, a system of polynomial equations is needed. This Julia package allows easy switching between a symbolic and numerical representation. For the symbolic object, construct a SymSys object as follows:

```julia
require("pnla.jl"); using PNLA;
eqs = :((x^2 + y^2 + z^2 - 1)/cos(0.3), x-y, 2/sqrt(3)*(x + y + y*z)^2*(-y))
var = :(x,y,z)
sys = SymSys(eqs, var)

eqst = "(x^2 + y^2 + z^2 - 1)/cos(0.3), x-y, 2/sqrt(3)*(x + y + y*z)*(-y)"
vart = "x,y,z"
syst = SymSys(eqst, vart)
```

Whereas division by polynomial other than a constant is not allowed, the user can easily use functions to compute coefficients. A factor based construction is allowed and simplifies designing systems based upon intersection and union of varieties using polynomials that vanish on the desired variety. The SymSys object holds in addition to the previously defined equations and variables, a Function object. It is recommended to use this function through the evalsys function in the PNLA package:

```julia
n = 3
A = rand(n,10)
evalsys(sys, A)
v1 = rand(3); v2 = rand(3); v3 = rand(3);
evalsys(sys, {v1,v2,v3})
evalsys(sys, v1, v2, v3)
```

Some special functions have been supplied through multiple dispatch. For example, evaluating a univariate polynomial system, and the specialized method can be accessed as follows:

```julia
methods(evalsys)
uni = :((x-1)*(x-1.5)*(x-3)*(x-4) + 10*x)
sys = SymSys(uni, :(x))
evalsys(sys, 2.0, 3.0, 4.0)
```

## Numerical representation

In order to start working with the numerical tools in the PNLA package, we need to construct a PolySys object out of the SymSys object. This can be achieved using the sym2poly function. In contrast to a SymSys object, which contains expression, the PolySys object has two fields:

1. Vector of coefficients for the Vandermonde monomials
2. The exponents of the Vandermonde monomials

The user is free to apply a monomial order, which will always be relative to the implicit order established by the variable expression found in the SymSys object. By default, the graded reverse lexicographic order is used:

```julia
polysys = sym2poly(sys)
polysys = sym2poly(sys, tdeg("x,y,z"))
polysys = sym2poly(sys, tdeg(:(x,y,z))

# Predefined monomial orders: lexicographic and graded lexicographic
polysys = sym2poly(sys, plex("y,x,z"))
polysys = sym2poly(sys, grlex("z,x,y"))

# Using a custom monomial order weighing matrix
polysys = sym2poly(sys, MonomialOrder([ 1 1 1; -1 0 0; 0 0 -1]))
```

## Symbolic representation

The reverse process requires a user-supplied list of variables. The method poly2sym expresses the PolySys object as a sum of monomials

```julia
symsys = poly2sym(polysys, "x,y,z")
symsys = poly2sym(polysys, :(x,y,z))
```

The user may supply a monomial order or opt for anonymous variables

```julia
symsys = poly2sym(polysys, "x,y,z", tdeg("y,x,z"))
symsys = poly2sym(polysys)
```

For better stability, the nested Horner scheme can be used. This will alter the construction of the equations expression. The same calls as with the poly2sym method can be made.

```julia
nestedsys = poly2horner(polysys, "x,y,z")
nestedsys = poly2horner(polysys)
```

## Basic operations

The prevalence of sparse polynomial expressions in many variables has led to an implementation using sparse matrix representations and parallel computing where deemed useful. Basic operations operate first and foremost on PolySys objects. For example

```julia
# Retrieve degree of each equation
deqeqs(polysys)
# Retrieve maximum degree of system
degmax(polysys)
# Bezout number (fails if overdetermined or underdetermined
bezout(polysys)
```
