pnla
====

Polynomial Numerical Linear Algebra package for Julia

# Introduction

This pakcga is in development. For a great PNLA package in a MATLAB/Oktave environment, check out Kim Batselier's code project [here](https://github.com/kbatseli/PNLA_MATLAB_OCTAVE).

# Constructing systems

In order to start working, a system of polynomial equations is needed. This Julia package allows easy switching between a symbolic and numerical representation. For the symbolic object, construct a SymSys object as follows:

```
include("pnla.jl"); using PNLA;
eqs = :(x^2 + y^2 + z^2 -1, x-y, 2*(x + y + y*z)*(-y))
var = :(x,y,z)
sys = SymSys(eqs, var)
```

The SymSys object holds in addition to the previously defined equations and variables, a Function object. It is recommended to use this function through the evalsys function in the PNLA package:

```
n = 3
A = rand(n,10)
evalsys(sys, A)
v1 = rand(3); v2 = rand(3); v3 = rand(3);
evalsys(sys, {v1,v2,v3})
evalsys(sys, v1, v2, v3)
```

Some special functions have been supplied through multiple dispatch. For example, evaluating a univariate polynomial system, and the specialized method can be accessed as follows:

```
methods(evalsys)
uni = :((x-1)*(x-1.5)*(x-3)*(x-4) + 10*x)
sys = SymSys(uni, :(x))
evalsys(sys, 2.0, 3.0, 4.0)
```
