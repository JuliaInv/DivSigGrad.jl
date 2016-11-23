using PyPlot
using jInv.Mesh
using jInv.Vis
using DivSigGrad
using KrylovMethods
using Base.Test


#
function fictousSourceTest2D(u,sig,rhs)
  refineMesh(M::RegularMesh) = (getRegularMesh(M.domain, M.n*2))

  M  = getRegularMesh([0 1 0 1],[3,4])
  N  = 6
  err = zeros(N,2)
  # k=3
  uk = 0.0; sk= 0.0; rk = 0.0;k=1;ut=0
   for k=1:N
      M  = refineMesh(M)
      xc = getCellCenteredGrid(M)
  	  xn = getNodalGrid(M)

  		uk = u(xn[:,1],xn[:,2])
      sk = sig(xc[:,1],xc[:,2])
      rk = rhs(xn[:,1],xn[:,2])

      A  = getDivSigGradMatrix(vec(sk),M)
  		V  = getVolume(M)
  		# G = getNodalGradientMatrix(M)
  		ut = A\(-vec(V[1,1]*rk))
  		# At = G'*G
  		# ut = cg(At,-vec(rk),tol=1e-15,maxIter=2000)[1]
  	ut = ut - (sum(ut)/prod(size(ut)))

  	res = norm(A*ut + rk,Inf);

      err[k,1] = V[1,1]*dot((ut-uk),A*(ut-uk))
      err[k,2] = norm(ut-uk,Inf)

      @printf "k=%d, n=[%d,%d], l2_err=%1.3e, factor=%1.3f linf_err=%1.3e\n" k M.n[1] M.n[2] err[k,1] err[max(k-1,1),1]/err[k,1] err[k,2]
   end

  @test countnz(diff(log(err[:,1])).<-1.8) > 4
end

# Constant sigma
fictousSourceTest2D((x,y)->cos(pi*x).*cos(pi*y), (x,y) -> ones(size(x)),
      (x,y) -> (- pi.^2.*cos(pi*x).*cos(pi*y) - pi.^2.*cos(pi*x).*cos(pi*y)))

# # VARYING SIGMA TEST
fictousSourceTest2D((x,y)->cos(pi*x).*cos(pi*y), (x,y) -> x.^2 + y.^2 + 1,
      (x,y) -> (- 2*pi^2*cos(pi*x).*cos(pi*y).*(x.^2 + y.^2 + 1) - 2*pi*x.*cos(pi*y).*sin(pi*x) - 2*pi.*y.*cos(pi*x).*sin(pi*y)))

println("\t== passed ! ==")
