import jInv.ForwardShare.getData
function getData(m::Vector{Float64},pFor::DivSigGradParam,doClear::Bool=false)
	A = getDivSigGradMatrix(m,pFor.Mesh)

	pFor.Ainv.doClear=1
	pFor.Fields,pFor.Ainv = solveLinearSystem(A,pFor.Sources,pFor.Ainv)		
	pFor.Ainv.doClear=0
	
	D = pFor.Receivers'*pFor.Fields
	return D, pFor
end