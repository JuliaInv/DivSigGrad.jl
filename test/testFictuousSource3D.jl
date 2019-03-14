#using PyPlot
using jInv.Mesh
using DivSigGrad
using KrylovMethods
using Test
using jInv.Utils

function fictousSourceTest3D(M,u,sig,rhs,expOrder=1.9)
  refineMesh(M::RegularMesh)  = (getRegularMesh(M.domain, M.n*2))
  refineMesh(M::TensorMesh3D) = (h1=rand(2*M.n[1]);h2=rand(2*M.n[2]);h3=rand(2*M.n[3]); getTensorMesh3D(h1/sum(h1),h2/sum(h2),h3/sum(h3)))
  N  = 5

  err = zeros(N,2)
  @printf "\t----Fictous Source Test DivSigGrad (%s)---\n\tk\tn\t\tl2_err\t\tfactor\tlinf_err\tfactor\n" typeof(M)
  uk = 0.0; sk= 0.0; rk = 0.0;k=1;ut=0
  for k=1:N
      M  = refineMesh(M)
      xc = getCellCenteredGrid(M)
      xn = getNodalGrid(M)
      
      # discretize solution, source, conductivity
      uk = u(xn[:,1],xn[:,2],xn[:,3])
      sk = sig(xc[:,1],xc[:,2],xc[:,3])
      rk = rhs(xn[:,1],xn[:,2],xn[:,3])

      # build PDE operator
      A  = getDivSigGradMatrix(vec(sk),M)
      V  = getVolume(M)
      An = getNodalAverageMatrix(M)
      W  = diagm(0 => vec(sum(An'*V,dims=2)))
      
      # solve PDE
      ut  = A\(-W*vec(rk))
      ut .-= mean(ut)

      # compute error
      err[k,1] = sqrt(dot((ut-uk),A*(ut-uk)))
      err[k,2] = norm(ut-uk,Inf)

      @printf "\t%d, n=[%-2d,%-2d,%-2d]\t\t%1.3e\t%1.3f\t%1.3e\t%1.3f\n" k M.n[1] M.n[2] M.n[2] err[k,1] err[max(k-1,1),1]/err[k,1] err[k,2] err[max(k-1,1),2]/err[k,2] 
	  if count(!iszero,diff(log2.(err[1:k,1])).<-expOrder) >= N-3
		break
	end
   end
  @test count(!iszero,diff(log2.(err[:,1])).<-expOrder) >= N-3
end

Mreg = getRegularMesh([0 1 0 1 0 1],[1,1,2])
h = rand(1); h /=sum(h)
Mten = getTensorMesh3D(h,h,h)
fictousSourceTest3D(Mreg,
							(x,y,z)-> cos.(pi*x) .* cos.(3*pi*y) .* cos.(2*pi*z),
							(x,y,z)-> tanh.(10*exp.(- 20*(x .- 1/2).^2 - 50*(y .- 1/2).^2 - 40*(z .- 1/2).^2)) .+ 1,
							(x,y,z)-> (- 14 .* pi.^2 .* cos.(pi .* x) .* cos.(3 .* pi .* y) .* cos.(2 .* pi .* z) .* 
							             (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 40 .* 
										 (z .- 1/2).^2)) .+ 1) - 10 .* pi .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* 
										 (y .- 1/2).^2 - 40 .* (z .- 1/2).^2) .* cos.(3 .* pi .* y) .* cos.(2 .* pi .* z) .*
										  sin.(pi .* x) .* (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 -
										   40 .* (z .- 1/2).^2)).^2 .- 1) .* (40 .* x .- 20) - 30 .* pi .* 
										   exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 40 .* (z .- 1/2).^2) .* 
										   cos.(pi .* x) .* cos.(2 .* pi .* z) .* sin.(3 .* pi .* y) .*
										    (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 40 .* 
											(z .- 1/2).^2)).^2 .- 1) .* (100 .* y .- 50) - 20 .* pi .* exp.(- 20 .* (x .- 1/2).^2 - 
											50 .* (y .- 1/2).^2 - 40 .* (z .- 1/2).^2) .* cos.(pi .* x) .* cos.(3 .* pi .* y) .* 
											sin.(2 .* pi .* z) .* (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 
											40 .* (z .- 1/2).^2)).^2 .- 1) .* (80 .* z .- 40)))

#
fictousSourceTest3D(Mten,
 									(x,y,z)-> cos.(pi*x) .* cos.(3*pi*y) .* cos.(2*pi*z),
 									(x,y,z)-> tanh.(10*exp.(- 20*(x .- 1/2).^2 - 50*(y .- 1/2).^2 - 40*(z .- 1/2).^2)) .+ 1,
 									(x,y,z)-> (- 14 .* pi.^2 .* cos.(pi .* x) .* cos.(3 .* pi .* y) .* cos.(2 .* pi .* z) .*
									           (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 40 .* 
											   (z .- 1/2).^2)) .+ 1) - 10 .* pi .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* 
											   (y .- 1/2).^2 - 40 .* (z .- 1/2).^2) .* cos.(3 .* pi .* y) .* cos.(2 .* pi .* z) .* 
											   sin.(pi .* x) .* (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 40 .* 
											   (z .- 1/2).^2)).^2 .- 1) .* (40 .* x .- 20) - 30 .* pi .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* 
											   (y .- 1/2).^2 - 40 .* (z .- 1/2).^2) .* cos.(pi .* x) .* cos.(2 .* pi .* z) .* 
											   sin.(3 .* pi .* y) .* (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 
											   40 .* (z .- 1/2).^2)).^2 .- 1) .* (100 .* y .- 50) - 20 .* pi .* exp.(- 20 .* (x .- 1/2).^2 -
											   50 .* (y .- 1/2).^2 - 40 .* (z .- 1/2).^2) .* cos.(pi .* x) .* cos.(3 .* pi .* y) .* 
											   sin.(2 .* pi .* z) .* (tanh.(10 .* exp.(- 20 .* (x .- 1/2).^2 - 50 .* (y .- 1/2).^2 - 40 .* 
											   (z .- 1/2).^2)).^2 .- 1) .* (80 .* z .- 40)), 1.5)

println("\t== passed ! ==")
