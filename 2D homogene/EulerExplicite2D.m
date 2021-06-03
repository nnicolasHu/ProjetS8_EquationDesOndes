function [t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny)
  dt=(EDP.T-EDP.t0)/Nt;
  dx=(EDP.b-EDP.a)/Nx;
  dy=(EDP.b-EDP.a)/Ny;
  
  t=EDP.t0 + dt*[0:Nt];
  x=EDP.a + dx*[0:Nx];
  y=EDP.c + dy*[0:Ny];
  
  %initialisation
  u=zeros((Ny+1)*(Nx+1),Nt+1); %le carré en espace est "compresse" en un vecteur colone
  u(:,1) = EvalFun2DVec(EDP.u0,x,y,Nx+1,Ny+1);
  Lap2D = Lap2DAssembling(Nx+1,Ny+1,dx,dy);
  u(:,2) = u(:,1) + dt*EvalFun2DVec(EDP.u1,x,y,Nx+1,Ny+1) + (dt*EDP.c)^2 /2 * Lap2D*u(:,1);
  
  
endfunction