
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



@testset "DivSigGrad" begin
@testset "getData" begin
	include("testGetData.jl")
end
@testset "DivSigGradMatrix" begin
	include("testDivSigGradMatrix.jl")
end
@testset "rectangular" begin
	include("testDivSigGrad.jl")
end
@testset "tensor" begin
	include("testDivSigGradTensor.jl")
end
@testset "fictous source (2D)" begin
	include("testFictuousSource2D.jl")
end
@testset "fictous source (3D)" begin
	include("testFictuousSource3D.jl")
end
end