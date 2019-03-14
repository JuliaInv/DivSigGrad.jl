using DivSigGrad
using Test

n1 = 8; n2 = 12; n3 = 9
x0    = [0.05, 0.1, 0.2]
domain = [x0[1], 1.2, x0[2], 1.4, x0[3], 2.6]

Mr = getRegularMesh(domain,[n1,n2,n3])

h1 = Mr.h[1].*ones(n1)
h2 = Mr.h[2].*ones(n2)
h3 = Mr.h[3].*ones(n3)
Mt = getTensorMesh3D(h1,h2,h3,x0)

m  = rand(n1,n2,n3)
Ar = getDivSigGradMatrix(vec(m),Mr)
At = getDivSigGradMatrix(vec(m),Mt)

@test norm(Ar-At,1)/norm(At,1) < 1e-13