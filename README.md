[![Build Status](https://travis-ci.org/JuliaInv/DivSigGrad.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/DivSigGrad.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/DivSigGrad.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/DivSigGrad.jl?branch=master)

# DivSigGrad.jl
Julia Package for Inverse Conducitivy Problems

# Requirements
This package is an add-on for [`jInv`](https://github.com/JuliaInv/jInv.jl), which needs to be installed. 

To accelerate the PDE solves it is also recommended to install [`ParSpMatVec`](https://github.com/lruthotto/ParSpMatVec) when using iterative solvers or [`MUMPS`](https://github.com/JuliaSparse/MUMPS.jl) as a direct solver. If these modules are available, they are used by default.

# Installation

In julia type:

``` 
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/DivSigGrad.jl","DivSigGrad")
Pkg.test("DivSigGrad")
```
