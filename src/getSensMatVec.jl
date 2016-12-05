import jInv.ForwardShare.getSensMatVec

function getSensMatVec(x::Vector{Float64},sig::Vector{Float64},pFor::DivSigGradParam)
      
    A  = getDivSigGradMatrix(sig,pFor.Mesh)
    G  = getNodalGradientMatrix(pFor.Mesh)
    Ae = getEdgeAverageMatrix(pFor.Mesh) 
    V  = getVolume(pFor.Mesh)
    Z  = G'*(sdiag(Ae'*V*x)*(G*pFor.Fields))
    Z, = solveLinearSystem(A,Z,pFor.Ainv)
    Jv = - getLinearMeasurements(pFor.Receivers,Z)
    return vec(Jv)
end