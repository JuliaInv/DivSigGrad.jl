import jInv.ForwardShare.getData
export getLinearMeasurements

"""
function DivSigGrad.getLinearMeasurements(Receivers,Fields)

get data by taking linear measurements from fields, e.g.,

	Receivers' * Fields

Depending on the experiment, Receivers can be a sparse matrix
(then same receivers are used for all fields) or an array of sparse
matrices (when receivers for each field are different)

Input:

	Receivers
	Fields

"""
function getLinearMeasurements(Receivers::SparseMatrixCSC,Fields::Array)
	return Receivers'*Fields
end

function getLinearMeasurements(Receivers::Array{SparseMatrixCSC},Fields::Array)
	nrec = length(Receivers)
	if size(Fields,2) != nrec
		error("There must be exactly one receiver matrix per field.")
	end

	D = zeros(size(Receivers[1],2),nrec)
	for k=1:nrec
		D[:,k] = Receivers[k]'*Fields[:,k]
	end
	return D
end


function getData(sig::Vector{Float64},pFor::DivSigGradParam,doClear::Bool=true)
	
	A = getDivSigGradMatrix(sig,pFor,doClear)

	pFor.Ainv.doClear=1
	pFor.Fields,pFor.Ainv = solveLinearSystem(A,pFor.Sources,pFor.Ainv)
	pFor.Ainv.doClear=0

	D = getLinearMeasurements(pFor.Receivers,pFor.Fields)
	return D, pFor
end
