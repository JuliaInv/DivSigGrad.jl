export getDivSigGradMatrix



"""
function DivSigGrad.getDivSigGradMatrix
	
builds finite volume discretization of diffusion operator

	A = DIV* ( sigma(x) GRAD u(x)),
	
where sigma is assumed to be a scalar quantity. Here, we use
Neuman boundary conditions and fix the value at u(0,0) to obtain
a regular matrix.

Inputs:

	sig::Vector{Float64}  - conductivities (cell-centered)
	M::AbstractTensorMesh - mesh from jInv.Mesh
	
Output:
	
	A  - PDE operator (sparse matrix)

"""
function getDivSigGradMatrix(sig::Vector{Float64},M::AbstractTensorMesh)
 	G       = getNodalGradientMatrix(M)
	Ae      = getEdgeAverageMatrix(M) 
	V       = getVolume(M)
    A       = (G'*sdiag(Ae'*(V*vec(sig))))*G
 	A[1,1] += 1
	return A
end


function getDivSigGradMatrix(sig::Vector{Float64},pFor::DivSigGradParam,doClear::Bool=true)
	if doClear || isempty(pFor.A) || size(pFor.A,1)!=prod(pFor.Mesh.n .+ 1)
		pFor.A = getDivSigGradMatrix(sig,pFor.Mesh)
	end
	return pFor.A
end