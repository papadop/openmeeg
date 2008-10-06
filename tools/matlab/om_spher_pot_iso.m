function u=om_spher_pot_iso(x,q,p,sigmas,rk) ;
% Calculates the electrical potential at point x on the surface of a layered
% conductive sphere induced by a current dipole of strength q at point
% p. The radii of the spheres are given by a vector rk in ascending order,
% the corresponding conductivities are in sigmas. SI units.
% The norm of x must be equal to the last rk.
%
% Usage: u=SpherPot(x,q,p,sigmas,rk) ;
%
% If called without arguments, makes a demo
%
% All vectors (x,q,p) should be in 3D.
%
% Based on formulas (1I),(2I) in Zhi Zhang: "A fast method to compute
% surface potential generated by dipoles within multilayer anisotropic
% spheres", in Phys. Med. Biol. 40 (1995), pages 335-349.
%
%
  % PROBLEM: Getting imaginary numbers solved ad hoc by real()

[n,m] = size(x);
if m~=3,
  error('x should be a set of row 3D vectors.')  
end

 q=q(:) ; p=p(:) ; % make vectors
q = q';p = p';
matp = ones(n,1)*p;
matq = ones(n,1)*q;

N=length(sigmas) ; % number of levels

% calculate angles 
alpha=om_vangle(matp,matq) ; cosgamma=cos(om_vangle(matp,x)) ;
beta=om_vangle(cross(matp,matq),cross(matp,x)) ;
relerr=1e-20 ; % keep adding terms as longer as the term to be added is greater

n=1 ; sum=0 ; 

while 1,
  term=(2*n+1)/n*(norm(p)/rk(N))^(n-1) ;
  % calculate matrix M
  M=eye(2)/(2*n+1)^(N-1) ;
  for k=1:N-1,
    rsigma=sigmas(k)/sigmas(k+1) ;
    rerk=(rk(N)/rk(k))^(2*n+1) ;
    M=M*[ n+(n+1)*rsigma (n+1)*(rsigma-1)*rerk ; ...
      n*(rsigma-1)/rerk (n+1)+n*rsigma ] ;
  end ;
  f=n/(n*M(2,2)+(1+n)*M(2,1)) ;
  P=legendre(n,cosgamma) ;
  % Attention! Matlab's definition of Legendre includes (-1)^m, unlike
  % the definition supposed for the following formula
  term=term*f*(n*cos(alpha).*P(1,:)'-cos(beta).*sin(alpha).*P(2,:)') ;
%   disp(['Adding term n=' num2str(n) ' value=' num2str(term) ]) ;
  sum=sum+term ;
  % finish if contribution is small or if something goes wrong
   if ((abs(term)<abs(sum)*relerr) & n>20)  | n>200,
    break ; 
  end ; 
  n=n+1 ;
end ;
u=sum*norm(q)/(4*pi*sigmas(N)*rk(N)^2) ;
u=real(u) ; % somewhere we were getting imaginary numbers, don't know why
% (with very small imaginary component)
