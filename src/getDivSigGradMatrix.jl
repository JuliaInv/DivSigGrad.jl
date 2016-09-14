export getDivSigGradMatrix



function getDivSigGradMatrix(m::Vector{Float64},Mesh::TensorMesh3D)
 	G       = getNodalGradientMatrix(Mesh)
	Ae      = getEdgeAverageMatrix(Mesh) 
	V       = getVolume(Mesh)
    A       = (G'*sdiag(Ae'*V*vec(m)))*G
 	A[1,1] += 1/Mesh.h1[1] 
	return A
end

function getDivSigGradMatrix(m::Vector{Float64},Mesh::RegularMesh)
 	G       = getNodalGradientMatrix(Mesh)
	Ae      = getEdgeAverageMatrix(Mesh) 
	V       = getVolume(Mesh)
    A       = (G'*sdiag(Ae'*V*vec(m)))*G
 	A[1,1] += 1/Mesh.h[1] 
	return A
end