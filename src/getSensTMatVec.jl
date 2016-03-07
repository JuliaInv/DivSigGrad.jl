import jInv.ForwardShare.getSensTMatVec

function getSensTMatVec(x::Vector{Float64},m::Vector{Float64},pFor::DivSigGradParam)
    
    A   = getDivSigGradMatrix(m,pFor.Mesh)
    U   = pFor.Fields
    G   = getNodalGradientMatrix(pFor.Mesh)
    Ae  = getEdgeAverageMatrix(pFor.Mesh) 
    x  = reshape(x,size(pFor.Receivers,2),size(pFor.Sources,2))
    
    Zt,  = solveLinearSystem(A,pFor.Receivers*x,pFor.Ainv)  
    JTv  = - sum(Ae*((G*U).*(G*Zt)),2)
    return vec(JTv)
end