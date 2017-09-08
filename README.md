[![Build Status](https://travis-ci.org/JuliaInv/DivSigGrad.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/DivSigGrad.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/DivSigGrad.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/DivSigGrad.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/rkal2ya6u0vd0vl0?svg=true)](https://ci.appveyor.com/project/lruthotto/divsiggrad-jl)


# DivSigGrad.jl

jInv Extension for solving Inverse Conductivity Problems. It supports various mesh types (regular meshes, stretched tensor meshes, OcTree meshes), different
PDE solvers (direct and iterative) and allows to easily parallelize inversions when multiple sources are available.

# Requirements

This package is intended to use with julia versions 0.6.x.

This package is an add-on for [`jInv`](https://github.com/JuliaInv/jInv.jl), which needs to be installed.

To accelerate the PDE solves, it is also recommended to install [`ParSpMatVec`](https://github.com/lruthotto/ParSpMatVec) when using iterative solvers or [`MUMPS`](https://github.com/JuliaSparse/MUMPS.jl) as a direct solver. If these modules are available, they are used by default.

# Installation

To use this package, start julia and type:

```
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/DivSigGrad.jl","DivSigGrad")
Pkg.test("DivSigGrad")
```

# Examples

A 3D inversion example can be found [`here`](https://github.com/JuliaInv/jInv.jl/blob/master/examples/exDCResistivity.ipynb).
