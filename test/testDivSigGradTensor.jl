using DivSigGrad
using jInv.Mesh
using jInv.LinearSolvers
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

Q = spzeros(Float64,prod(M.n+1),6)
cnt = 1
for i=5:5:20
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



fields = [0.0] 
@everywhere PCGsolver(A,b,M;kwargs...) = cg(A,b,M=M;kwargs...)
Apcg         = getIterativeSolver(PCGsolver)
Apcg.tol=1e-10
Abpcg      = getBlockIterativeSolver(KrylovMethods.blockCG,tol=1e-5)
Abpcg.out=0
Abpcg.maxIter = 50000
Abpcg.tol  = 1e-10

pFors     = [];
push!(pFors,DivSigGradParam(M,Q,P,fields,Apcg))
push!(pFors,DivSigGradParam(M,Q,P,fields,Abpcg))

# different receivers for each source
Rec = Array{SparseMatrixCSC}(size(Q,2))
for k=1:size(Q,2)
	Rec[k] = P
end
push!(pFors,DivSigGradParam(M,Q,Rec,fields,Apcg))
push!(pFors,DivSigGradParam(M,Q,Rec,fields,Abpcg))

if LinearSolvers.hasMUMPS
	Amumps    = getMUMPSsolver()
	push!(pFors,DivSigGradParam(M,Q,P,fields,Amumps))
end

# Forward problem
m = ones(n1,n2,n3);
m[round(Int64,n1/3):round(Int64,n1/2),
  round(Int64,n1/3):round(Int64,n1/2),
  end-4:end-2] = 2

(D,pFors[1]) = getData(m[:],pFors[1]);
D0, = getData(m[:]*0+1.0,pFors[1]);

for k=2:length(pFors)
	
	(Db,pFors[k]) = getData(m[:],pFors[k]);
	Db0, = getData(m[:]*0+1.0,pFors[k]);
	@test norm(D-Db)/norm(D) < 1e-1
	@test norm(D0-Db0)/norm(D0) < 1e-1
end

# Derivative check
for k=1:length(pFors)
	println("\t--- derivative for solver $(typeof(pFors[k].Ainv)) ---")
	dm = randn(size(m))*1e-1
	Jdm = getSensMatVec(dm[:],m[:],pFors[k])
	alpha = 1.0;
	err = zeros(6,2)
	for i=1:size(err,1)
		(D1,pFors[k]) = getData(m[:]+alpha*dm[:],pFors[k]);
   		err[i,1] = norm(D1[:]-D[:])
   		err[i,2] = norm(D1[:]-D[:]-alpha*Jdm)
		@printf "\talpha=%1.2e\t\tE0=%1.2e\t\tE1=%1.2e\n" alpha err[i,1] err[i,2]
		alpha = alpha/2
	end
	@test length(find(2+diff(log2(err[:,2])).<0.2))>=3
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
