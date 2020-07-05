using SparseArrays
using jInv.Mesh
using DivSigGrad
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.Utils
using Test
using KrylovMethods

n1    = 20; n2 = 18; n3 = 8;
ns    = 100;
M     = getRegularMesh([0,1,0,1,0,0.5],[n1,n2,n3])

Q = spzeros(Float64,prod(M.n .+ 1),6)
cnt = 1
for i=5:2:15
	global cnt
	q     = zeros(n1+1,n2+1,n3+1);
	q[i,2,end]     = 1.0/prod(M.h)
	q[i,end-1,end] = -1.0/prod(M.h)
	Q[:,cnt]       .= sparse(q[:])
	cnt+=1
end

P     = spzeros((n1+1)*(n2+1)*(n3+1),(n1)*(n2));
cnt   = 1
for i=1:n1
	for j=1:n2
		global cnt
		p        = zeros(n1+1,n2+1,n3+1);
      p[i,j,end] = 1.0/prod(M.h); p[i+1,j,end] = -1.0/prod(M.h)
		P[:,cnt] .= sparse(p[:])
		cnt += 1
	end
end

fields = [0.0]
@everywhere PCGsolver(A,b;M=M,tol=1e-8,maxIter=1000,out=-1) = KrylovMethods.cg(A,b;M=M,tol=tol,maxIter=maxIter,out=out)
Apcg         = getIterativeSolver(PCGsolver)
Apcg.maxIter=1000
@everywhere IterMethod(A,B;M=M,X=X,tol=1e-8,maxIter=50000,out=-1) = KrylovMethods.blockCG(A,B;M=M,X=X,tol=tol,maxIter=maxIter,out=out)
Abpcg      = getBlockIterativeSolver(IterMethod);
Abpcg.out=0
Abpcg.maxIter = 50000
Abpcg.tol  = 1e-8


Ppcg      = getDivSigGradParam(M,Q,P,Ainv=Apcg)
Pbpcg     = getDivSigGradParam(M,Q,P,Ainv=Abpcg)

Ajulia    = getJuliaSolver()
Pjulia    = getDivSigGradParam(M,Q,P,Ainv=Ajulia)


# Forward problem
m = ones(n1,n2,n3);
m[round(Int64,n1/3):round(Int64,n1/2),
  round(Int64,n1/3):round(Int64,n1/2),
  end-4:end-2] .= 2

println("use PCG")
@time (D,Ppcg) = getData(m[:],Ppcg);
D0, = getData(m[:]*0 .+ 1.0,Ppcg);

println("use BlockPCG")
@time (Db,Pbpcg) = getData(m[:],Pbpcg);
Db0, = getData(m[:]*0 .+ 1.0,Pbpcg);
@test norm(D-Db)/norm(D) < 1e-1
@test norm(D0-Db0)/norm(D0) < 1e-1


println("use Julia Solver")
@time (Dm,Pjulia) = getData(m[:],Pjulia);
D0m, = getData(m[:]*0 .+ 1.0,Pjulia);
@test norm(D-Dm)/norm(Dm) < 1e-1
@test norm(D0-D0m)/norm(D0m) < 1e-1


# Derivative check
println("\t--- derivative for PCG ---")
(Jn,Jm) = getSensMatSize(Ppcg)
pass, = checkDerivative(vec(m),Ppcg)
@test pass

println("\t--- derivative for BlockPCG ---")
pass, = checkDerivative(vec(m),Pbpcg)
@test pass

# Derivative check

println("\t--- derivative Julia Solver ---")
pass, = checkDerivative(vec(m),Pjulia)
@test pass


println("\t--- adjoint test PCG ---")
dm = randn(Jm)*1e-1
Jdm       = getSensMatVec(dm[:],m[:],Ppcg)
v         = randn(Jn)
t1        = dot(v,Jdm)
JTv       = getSensTMatVec(v[:],m[:],Ppcg)
t2        = dot(JTv,dm[:])
@test (abs(t1-t2)/t1<=5e-2)

println("\t--- adjoint test BlockPCG ---")
Jdm = getSensMatVec(dm[:],m[:],Pbpcg)
v         = randn(Jn)
t1        = dot(v,Jdm)
JTv       = getSensTMatVec(v[:],m[:],Pbpcg)
t2        = dot(JTv,dm[:])
@test (abs(t1-t2)/t1<=5e-2)


println("\t--- adjoint test Julia Solver ---")
Jdm        = getSensMatVec(dm[:],m[:],Pjulia)
v         = randn(Jn)
t1        = dot(v,Jdm)
JTv       = getSensTMatVec(v[:],m[:],Pjulia)
t2        = dot(JTv,dm[:])
@test (abs(t1-t2)/t1<=5e-2)
