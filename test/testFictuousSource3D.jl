#using PyPlot
using jInv.Mesh
using jInv.Vis
using DivSigGrad
using KrylovMethods
using Base.Test
using jInv.Utils



function fictousSourceTest3D(M,u,sig,rhs,expOrder=2.5)
	refineMesh(M::RegularMesh)  = (getRegularMesh(M.domain, M.n*2))
  refineMesh(M::TensorMesh3D) = (h1=rand(2*M.n[1]);h2=rand(2*M.n[2]);h3=rand(2*M.n[3]); getTensorMesh3D(h1/sum(h1),h2/sum(h2),h3/sum(h3)))
  N  = 4

	err = zeros(N,2)
  # k=3
  uk = 0.0; sk= 0.0; rk = 0.0;k=1;ut=0
   for k=1:N
      M  = refineMesh(M)
      xc = getCellCenteredGrid(M)
  	  xn = getNodalGrid(M)

  		uk = u(xn[:,1],xn[:,2],xn[:,3])
      sk = sig(xc[:,1],xc[:,2],xc[:,3])
      rk = rhs(xn[:,1],xn[:,2],xn[:,3])

      A  = getDivSigGradMatrix(vec(sk),M)
  		V  = getVolume(M)
  		An = getNodalAverageMatrix(M)
			W  = diagm(vec(sum(An'*V,2)))
			ut = A\(-W*rk)

		  ut = ut - mean(ut)

  	res = norm(A*ut + rk,Inf);

      err[k,1] = sqrt(dot(W*(ut-uk),A*(ut-uk)))
      err[k,2] = norm(ut-uk,Inf)

      @printf "k=%d, n=[%d,%d,%d], l2_err=%1.3e, factor=%1.3f linf_err=%1.3e\n" k M.n[1] M.n[2] M.n[2] err[k,1] err[max(k-1,1),1]/err[k,1] err[k,2]
   end

	#  figure(1); clf()
	# println(diff(log2(err[:,1])))
	#  viewOrthoSlices2D(sk,M)
  @test countnz(diff(log2(err[:,1])).<-expOrder) >= N-2
end

Mreg = getRegularMesh([0 1 0 1 0 1],[1,1,2])
h = rand(1); h /=sum(h)
Mten = getTensorMesh3D(h,h,h)
fictousSourceTest3D(Mreg,
										(x,y,z)-> cos(pi*x).*cos(3*pi*y).*cos(2*pi*z),
										(x,y,z)-> tanh(10*exp(- 20*(x - 1/2).^2 - 50*(y - 1/2).^2 - 40*(z - 1/2).^2)) + 1,
										(x,y,z)-> (- 14.*pi.^2.*cos(pi.*x).*cos(3.*pi.*y).*cos(2.*pi.*z).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)) + 1) - 10.*pi.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2).*cos(3.*pi.*y).*cos(2.*pi.*z).*sin(pi.*x).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)).^2 - 1).*(40.*x - 20) - 30.*pi.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2).*cos(pi.*x).*cos(2.*pi.*z).*sin(3.*pi.*y).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)).^2 - 1).*(100.*y - 50) - 20.*pi.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2).*cos(pi.*x).*cos(3.*pi.*y).*sin(2.*pi.*z).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)).^2 - 1).*(80.*z - 40)))

#
fictousSourceTest3D(Mten,
 									(x,y,z)-> cos(pi*x).*cos(3*pi*y).*cos(2*pi*z),
 									(x,y,z)-> tanh(10*exp(- 20*(x - 1/2).^2 - 50*(y - 1/2).^2 - 40*(z - 1/2).^2)) + 1,
 									(x,y,z)-> (- 14.*pi.^2.*cos(pi.*x).*cos(3.*pi.*y).*cos(2.*pi.*z).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)) + 1) - 10.*pi.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2).*cos(3.*pi.*y).*cos(2.*pi.*z).*sin(pi.*x).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)).^2 - 1).*(40.*x - 20) - 30.*pi.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2).*cos(pi.*x).*cos(2.*pi.*z).*sin(3.*pi.*y).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)).^2 - 1).*(100.*y - 50) - 20.*pi.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2).*cos(pi.*x).*cos(3.*pi.*y).*sin(2.*pi.*z).*(tanh(10.*exp(- 20.*(x - 1/2).^2 - 50.*(y - 1/2).^2 - 40.*(z - 1/2).^2)).^2 - 1).*(80.*z - 40)),
									1.5)

println("\t== passed ! ==")
