function [t,x,y,u]=EulerExpliciteNaive(EDP,nt,nx,ny)
  Axy=Lap2DAssembling(nx,ny,hx,hy);
  N=(nx)*(ny);
  
  dt=(EDP.T-EDP.t0)/Nt;
  dx=(EDP.b-EDP.a)/Nx;
  dy=(EDP.d-EDP.c)/Ny;
  
  t=EDP.t0 + dt*[0:Nt];
  x=EDP.a + dx*[0:Nx];
  y=EDP.c + dy*[0:Ny];
  
  cte=(EDP.c*dt)**2;
  
  %initialisation
  V=zeros((Nx+1)*(Ny+1),Nt+1); %colonne=espace et ligne=temps
  
  for k=1:N
    [i,j]=bijRecF(k);
    V(k,1)=EDP.u0(i,j)'; %bord gauche (t=t0)
    V(k,2)=( spMatDiag(ones(1,nx*ny)) + (cte/2)*Lap1D(nx*ny))*EDP.u0(i,j)' + dt*EDP.u1(i,j)';
  endfor
  
  V(1,:)=EDP.ua(t); %bord haut 
  V(end,:)=EDP.ub(t); %bord bas 
  
  for i=2:N-1
    V(i,n+1) = cte*(((V(i+1,n)-2*V(i,n)+V(i-1,N))/(dx*dx))+((V(i+nx,n)-2*V(i,n)+V(i-nx,N))/(dy*dy)))+2*V(i,n)-V(i,n-1);  
  endfor
  
  
endfunction
