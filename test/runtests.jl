
using Test
using LinearAlgebra
using Distributed
using Printf
using Statistics


if nworkers()<2
	addprocs(2)
end

using jInv.Mesh
using jInv.ForwardShare
using jInv.Utils
using DivSigGrad
using jInv.LinearSolvers
using KrylovMethods


println("==== test getData.jl for rectangular mesh ===")
include("testGetData.jl")
println("==== test DivSigGradMatrix for tensor and regular mesh ===")
include("testDivSigGradMatrix.jl")
println("==== test DivSigGrad.jl for rectangular mesh ===")
include("testDivSigGrad.jl")
println("==== test DivSigGrad.jl for tensor mesh ===")
include("testDivSigGradTensor.jl")
println("==== fictous source test (2D) ===")
include("testFictuousSource2D.jl")
println("==== fictous source test (3D) ===")
include("testFictuousSource3D.jl")
