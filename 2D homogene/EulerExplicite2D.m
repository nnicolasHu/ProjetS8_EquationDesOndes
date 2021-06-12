function [t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny)
  dt=(EDP.T-EDP.t0)/Nt;
  dx=(EDP.b-EDP.a)/Nx;
  dy=(EDP.b-EDP.a)/Ny;
  
  t=EDP.t0 + dt*[0:Nt];
  x=EDP.a + dx*[0:Nx];
  y=EDP.a + dy*[0:Ny];
  
  %initialisation
  u=zeros((Ny+1)*(Nx+1),Nt+1); %le carré en espace est "compresse" en un vecteur colone
  u(:,1) = EvalFun2DVec(EDP.u0,x,y,Nx+1,Ny+1);
  Lap2D = Lap2DAssembling(Nx+1,Ny+1,dx,dy);
  u(:,2) = u(:,1) + dt*EvalFun2DVec(EDP.u1,x,y,Nx+1,Ny+1) + (dt^2 /2)*(EvalFun2DVec(EDP.f,x,y,Nx+1,Ny+1) + EDP.c^2 * Lap2D*u(:,1));
  
  %Bord
  BD = BordDomaine(Nx+1,Ny+1);
  for n=2:Nt
    u(:,n+1) = Lap2D*u(:,n) + 2*u(:,n) - u(:,n-1);
  endfor
  
  
endfunction