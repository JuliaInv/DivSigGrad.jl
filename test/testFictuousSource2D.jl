#using PyPlot
using jInv.Mesh
using jInv.Utils
using DivSigGrad
using KrylovMethods
using Test


#
function fictousSourceTest2D(u,sig,rhs)
  refineMesh(M::RegularMesh) = (getRegularMesh(M.domain, M.n*2))

  M  = getRegularMesh([0 1 0 1],[3,4])
  N  = 6
  err = zeros(N,2)
  
  @printf "----Fictous Source Test DivSigGrad 2D---\nk\tn\t\tl2_err\t\tfactor\tlinf_err\tfactor\n"
  for k=1:N
      M  = refineMesh(M)
      xc = getCellCenteredGrid(M)
      xn = getNodalGrid(M)

      uk = u(xn[:,1],xn[:,2])
      sk = sig(xc[:,1],xc[:,2])
      rk = rhs(xn[:,1],xn[:,2])

      A  = getDivSigGradMatrix(vec(sk),M)
      V  = getVolume(M)
      v  = Vector(diag(V))
      An = getNodalAverageMatrix(M)
      ut = A\( (An'*v).*vec(-rk))
      ut .-= mean(ut)

      err[k,1] = sqrt(dot((ut-uk),A*(ut-uk)))
      err[k,2] = norm(ut-uk,Inf)

      @printf "%d\t[%-3d,%-3d]\t%1.3e\t%1.3f\t%1.3e\t%1.3f\n" k M.n[1] M.n[2] err[k,1] err[max(k-1,1),1]/err[k,1] err[k,2] err[max(k-1,1),2]/err[k,2]
   end
  @test count(!iszero,diff(log2.(err[:,1])).<-1.8) > 4
end

# Constant sigma
fictousSourceTest2D((x,y)->cos.(pi*x).*cos.(pi*y), (x,y) -> ones(size(x)),
      (x,y) -> (- pi.^2 .* cos.(pi*x).*cos.(pi*y) - pi.^2 .* cos.(pi*x).*cos.(pi*y)))

# # VARYING SIGMA TEST
fictousSourceTest2D((x,y)->cos.(pi*x).*cos.(pi*y), (x,y) -> x.^2 + y.^2 .+ 1,
      (x,y) -> (- 2*pi^2*cos.(pi*x).*cos.(pi*y).*(x.^2 + y.^2 .+ 1) - 2*pi*x.*cos.(pi*y).*sin.(pi*x) - 2*pi.*y.*cos.(pi*x).*sin.(pi*y)))

println("\t== passed ! ==")
