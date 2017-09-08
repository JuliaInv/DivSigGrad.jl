using DivSigGrad
using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.Utils
using Base.Test
using KrylovMethods

# construct tensor mesh
n1 = 20; n2 = 18; n3 = 8;
h1    = 1+.1*rand(n1)
h2    = 1+.1*rand(n2)
h3    = 1+.1*rand(n3)
ns    = 100
M     = getTensorMesh3D(h1,h2,h3)

srci  =5:5:20; 
Q = spzeros(Float64,prod(M.n+1),length(srci))
cnt = 1
for i=srci
	q     = zeros(n1+1,n2+1,n3+1);
	q[i,5,end] = 1.0
	q[i,end-5,end] = -1.0
	Q[:,cnt] = sparse(q[:])
	cnt+=1
end

P     = spzeros((n1+1)*(n2+1)*(n3+1),(n1)*(n2));
cnt   = 1
for i=1:n1
	for j=1:n2
		p     = zeros(n1+1,n2+1,n3+1);
      p[i,j,end] = 1.0; p[i+1,j,end] = -1.0
		P[:,cnt] = sparse(p[:])
		cnt += 1
	end
end

# @everywhere PCGsolver(A,b,M;kwargs...) = cg(A,b,M=M;kwargs...)
Apcg         = getIterativeSolver(cg)
Apcg.tol=1e-10
Abpcg         = getBlockIterativeSolver(KrylovMethods.blockCG,tol=1e-10)
Abpcg.out     =0
Abpcg.maxIter = 50000
Abpcg.tol     = 1e-10

pFors     = [];
push!(pFors,getDivSigGradParam(M,Q,P,Ainv=Apcg))
push!(pFors,getDivSigGradParam(M,Q,P,Ainv=Abpcg))

# different receivers for each source
Rec = Array{SparseMatrixCSC}(size(Q,2))
for k=1:size(Q,2)
	Rec[k] = P
end
push!(pFors,getDivSigGradParam(M,Q,Rec,Ainv=Apcg))
push!(pFors,getDivSigGradParam(M,Q,Rec,Ainv=Abpcg))

if LinearSolvers.hasMUMPS
	Amumps    = getMUMPSsolver()
	push!(pFors,getDivSigGradParam(M,Q,P,Ainv=Amumps))
end

# Forward problem
m = ones(n1,n2,n3);
m[round(Int64,n1/3):round(Int64,n1/2),
  round(Int64,n1/3):round(Int64,n1/2),
  end-4:end-2] = 2

(D,pFors[1]) = getData(m[:],pFors[1]);
D0, = getData(m[:]*0+1.0,pFors[1]);

for k=2:length(pFors)
	println("\t--- getData for solver $(typeof(pFors[k].Ainv)) ---")

	(Db,pFors[k]) = getData(m[:],pFors[k]);
	Db0, = getData(m[:]*0+1.0,pFors[k]);
	@test norm(D-Db)/norm(D) < 1e-1
	@test norm(D0-Db0)/norm(D0) < 1e-1
end

# Derivative check
for k=1:length(pFors)
	println("\t--- derivative for solver $(typeof(pFors[k].Ainv)) ---")
	pass, = checkDerivative(vec(m),pFors[k])
	@test pass
end

for k=1:length(pFors)
	println("\t--- adjoint test for solver $(typeof(pFors[k].Ainv)) ---")
	dm = randn(size(m))*1e-1
	Jdm = getSensMatVec(dm[:],m[:],pFors[k])
	v         = randn(size(Jdm))
	t1        = dot(v,Jdm)
	JTv       = getSensTMatVec(v[:],m[:],pFors[k])
	t2        = dot(JTv,dm[:])
	@test (abs(t1-t2)/t1<=5e-2)
end
