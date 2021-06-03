function [t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny)
  dt=(EDP.T-EDP.t0)/Nt;
  dx=(EDP.b-EDP.a)/Nx;
  dy=(EDP.d-EDP.c)/Ny;
  
  t=EDP.t0 + dt*[0:Nt];
  x=EDP.a + dx*[0:Nx];
  y=EDP.c + dy*[0:Ny];
  
  %initialisation
  u=zeros(Ny+1,Nx+1,Nt+1); %1er coord = y(ligne), 2� coord = x(colonne), 3� coord = temps(hauteur)
  u(:,:,1)=calculEspace(EDP.u0,x,y);
  
  
endfunction