import jInv.ForwardShare.getSensTMatVec

function getSensTMatVec(x::Vector{Float64},m::Vector{Float64},pFor::DivSigGradParam)
    
    A   = getDivSigGradMatrix(m,pFor.Mesh)
    U   = pFor.Fields
    G   = getNodalGradientMatrix(pFor.Mesh)
	V   = getVolume(pFor.Mesh)
    Ae  = getEdgeAverageMatrix(pFor.Mesh) 
	rhs = getRHS(pFor.Receivers,pFor.Sources,x)
    Zt,  = solveLinearSystem(A,rhs,pFor.Ainv)  
    JTv  = - sum(V*Ae*((G*U).*(G*Zt)),2)
    return vec(JTv)
end

function getRHS(Receivers::SparseMatrixCSC,Sources,data)
	data = reshape(data,size(Receivers,2),size(Sources,2))
	return Receivers*data
end

function getRHS(Receivers::Array{SparseMatrixCSC},Sources,data)
	nrec = length(Receivers)
	if size(Sources,2) != nrec
		error("There must be exactly one receiver matrix per source.")
	end
    data = reshape(data,size(Receivers[1],2),size(Sources,2))

	rhs  = zeros(size(Receivers[1],1),size(Sources,2))
	for k=1:nrec
		rhs[:,k] = Receivers[k]*data[:,k]
	end
	return rhs
end