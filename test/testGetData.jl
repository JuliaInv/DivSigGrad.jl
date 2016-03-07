using jInv.ForwardShare
@everywhere begin
	using jInv.Mesh
	using DivSigGrad
	using jInv.LinearSolvers
	using KrylovMethods
end
if nworkers()==1
	warn("add more workers to test parallel code")
end
using Base.Test

# get mesh for conductivity
ns    = 100;
Minv  = getRegularMesh([0,1,0,1,0,0.5],[32,32,16])

# get mesh for forward problems
n1    = 32; n2 = 32; n3 = 16;
Mfor  = getRegularMesh([0,1,0,1,0,0.5],[n1,n2,n3])

# first source
Q1     = zeros(n1+1,n2+1,n3+1);
Q1[round(Int,n1/4),round(Int,n1/2),end] = 1.0
Q1[round(Int,3*n1/4),round(Int,n1/2),end] = -1.0
Q1     = sparse(Q1[:])

# second source
Q2     = zeros(n1+1,n2+1,n3+1);
Q2[round(Int,3*n1/4),round(Int,n1/2),end] = 1.0
Q2[round(Int,3*n1/4),round(Int,n1/2),end] = -1.0
Q2     = sparse(Q2[:])

P     = spzeros((n1+1)*(n2+1)*(n3+1),(n1)*(n2));
cnt   = 1
for i=1:n1
	for j=1:n2
		p     = zeros(n1+1,n2+1,n3+1);
      p[i,j] = 1.0; p[i+1,j] = -1.0
		P[:,cnt] = sparse(p[:])
		cnt += 1
	end
end

fields = [0.0] 
# Ainv  = getMUMPSsolver()
@everywhere PCGsolver(A,b,M;kwargs...) = cg(A,b,M=M;kwargs...)
Ainv         = getIterativeSolver(PCGsolver)

# distribute forward problems
PF = Array{RemoteRef{Channel{Any}}}(2)
PF[1] = @spawn DivSigGradParam(Mfor,Q1,P,fields,Ainv)
PF[2] = @spawn DivSigGradParam(Mfor,Q2,P,fields,Ainv)

# get mesh to mesh interpolation
M2M   = prepareMesh2Mesh(PF,Minv,false)
m = rand(prod(Minv.n))
@time begin
	Dobs1t, = getData(fetch(M2M[1])'*vec(m),fetch(PF[1]))
	Dobs2t, = getData(fetch(M2M[2])'*vec(m),fetch(PF[2]))
end

@time begin
	Dobsr,PF = getData(vec(m),PF,M2M)
	Dobs1 = fetch(Dobsr[1])
	Dobs2 = fetch(Dobsr[2])
end
@test_approx_eq Dobs1 Dobs1t
@test_approx_eq Dobs2 Dobs2t


M2M   = prepareMesh2Mesh(PF,Minv,true)
@time begin
	Dobsr,PF = getData(vec(m),PF,M2M)
	Dobs1 = fetch(Dobsr[1])
	Dobs2 = fetch(Dobsr[2])
end
@test_approx_eq Dobs1 Dobs1t
@test_approx_eq Dobs2 Dobs2t




