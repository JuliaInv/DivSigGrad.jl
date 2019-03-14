if nworkers()<2
	addprocs(2)
end

using SparseArrays
using jInv.Mesh
using jInv.ForwardShare
using jInv.Utils
using DivSigGrad
using jInv.LinearSolvers
using KrylovMethods
using Test


# get mesh for conductivity
ns    = 100;
Minv  = getRegularMesh([0,1,0,1,0,0.5],[30,30,10])

# get mesh for forward problems
n1    = 15; n2 = 15; n3 = 5;
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
		global cnt
		p     = zeros(n1+1,n2+1,n3+1);
        p[i,j,1] = 1.0; p[i+1,j,1] = -1.0
		P[:,cnt] = sparse(p[:])
		cnt += 1
	end
end

fields = [0.0]
Ainv         = getIterativeSolver(KrylovMethods.cg)

# distribute forward problems
PF    = Array{RemoteChannel}(undef,2)
PF[1] = initRemoteChannel(DivSigGradParam,workers()[1],Mfor,[Q1 Q2],[P P],fields,spzeros(0,0),Ainv)
PF[2] = initRemoteChannel(DivSigGradParam,workers()[2],Mfor,[Q2 Q2],[P P],fields,spzeros(0,0),Ainv)

@test getNumberOfData(PF[1])==getNumberOfData(PF[2])
@test getNumberOfData(PF)==2*getNumberOfData(PF[2])
# get mesh to mesh interpolation
M2M    = prepareMesh2Mesh(PF,Minv,false)
m      = rand(prod(Minv.n))
sig    = exp.(m)
sigloc = fetch(M2M[1])'*vec(sig)

@time begin
	Dobs1t, = getData(sigloc,fetch(PF[1]))
	Dobs2t, = getData(sigloc,fetch(PF[2]))
end

@time begin
	Dobsr,PF = getData(vec(sig),PF,M2M)
	Dobs1 = fetch(Dobsr[1])
	Dobs2 = fetch(Dobsr[2])
end
@test Dobs1 ≈ Dobs1t
@test Dobs2 ≈ Dobs2t


M2M   = prepareMesh2Mesh(PF,Minv,true)
@time begin
	Dobsr,PF = getData(vec(sig),PF,M2M)
	Dobs1 = fetch(Dobsr[1])
	Dobs2 = fetch(Dobsr[2])
end
@test Dobs1 ≈ Dobs1t
@test Dobs2 ≈ Dobs2t


# test for different receivers for different sources
P1 = P[:,1:round(Int64,n1*n2/2)]     # receivers for first source
P2 = P[:,round(Int64,n1*n2/2)+1:end-1] # receivers for second source

# Option 1: put them in same param
Ps = Array{SparseMatrixCSC}(undef,2)
Ps[1] = P1
Ps[2] = P2
PF1 = getDivSigGradParam(Mfor,[Q1 Q2],Ps,Ainv=Ainv)
Dobs1, = getData(sigloc,PF1)

# Option 2: use two params
PF2 = Array{RemoteChannel}(undef,2)
PF2[1] = initRemoteChannel(DivSigGradParam, workers()[1], Mfor,Q1,P1,fields,spzeros(0,0),Ainv)
PF2[2] = initRemoteChannel(DivSigGradParam, workers()[2], Mfor,Q2,P2,fields,spzeros(0,0),Ainv)
Dobs2, = getData(vec(sig),PF2,M2M)

for k=1:length(Ps)
	@test Dobs1[:,k] ≈ fetch(Dobs2[k])
end
