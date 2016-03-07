import jInv.ForwardShare.getSensMatVec

function getSensMatVec(x::Vector{Float64},m::Vector{Float64},pFor::DivSigGradParam)
      
    A  = getDivSigGradMatrix(m,pFor.Mesh)
    G  = getNodalGradientMatrix(pFor.Mesh)
    Ae = getEdgeAverageMatrix(pFor.Mesh) 
    
    Z  = G'*(sdiag(Ae'*x)*(G*pFor.Fields))
    Z, = solveLinearSystem(A,Z,pFor.Ainv)
    Jv = -pFor.Receivers'*Z
    return vec(Jv)
end