module DivSigGrad
export getDivSigGradMatrix,DivSigGradParam, getData, getSensMatVec, getSensTMatVec

using jInv.Mesh
using jInv.LinearSolvers
using jInv.Utils


export DivSigGradParam
import jInv.ForwardShare.ForwardProbType

type DivSigGradParam <: ForwardProbType
    Mesh::AbstractMesh
    Sources::SparseMatrixCSC
    Receivers::SparseMatrixCSC
    Fields::Array{Float64} 
    Ainv::AbstractSolver
end

include("getDivSigGradMatrix.jl")	
include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")

end