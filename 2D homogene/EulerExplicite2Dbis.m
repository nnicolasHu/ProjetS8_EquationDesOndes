function [t,x,y,u]=EulerExplicite2Dbis(EDP,Nt,Nx,Ny)
  dt=(EDP.T-EDP.t0)/Nt;
  dx=(EDP.b-EDP.a)/Nx;
  dy=(EDP.b-EDP.a)/Ny;
  
  t=EDP.t0 + dt*[0:Nt];
  x=EDP.a + dx*[0:Nx];
  y=EDP.a + dy*[0:Ny];
  
##  %Bord
  BD = BordDomaine(Nx+1,Ny+1);
  idxBD=unique([BD.Sud,BD.Nord,BD.Est,BD.Ouest]);
  idxBDc=setdiff(1:(Ny+1)*(Nx+1),idxBD);
  [I,J]=bijRecF(1:(Ny+1)*(Nx+1),Nx+1);
  X=x(I+1);
  Y=y(J+1);

  %initialisation en temps
  u=zeros((Ny+1)*(Nx+1),Nt+1); %le carré en espace est "compresse" en un vecteur colone
  u(:,1) = EvalFun2DVec(EDP.u0,x,y,Nx+1,Ny+1);
  Lap2D = Lap2DAssembling(Nx+1,Ny+1,dx,dy);
  u(:,2) = u(:,1) + dt*EvalFun2DVec(EDP.u1,x,y,Nx+1,Ny+1) + ( dt^2 /2)*(EvalSource2DVec(EDP.f,x,y,Nx+1,Ny+1,t(1)) + EDP.c^2*Lap2D*u(:,1));
  u(BD.Nord,2)=EDP.ubord(t(2),X(BD.Nord),Y(BD.Nord))';
  u(BD.Sud,2)=EDP.ubord(t(2),X(BD.Sud),Y(BD.Sud))';
  u(BD.Est,2)=EDP.ubord(t(2),X(BD.Est),Y(BD.Est))';
  u(BD.Ouest,2)=EDP.ubord(t(2),X(BD.Ouest),Y(BD.Ouest))';
 
  for n=2:Nt
    u(:,n+1) = ((EDP.c*dt)^2)*Lap2D*u(:,n) + 2*u(:,n) - u(:,n-1)+ (dt^2)*EvalSource2DVec(EDP.f,x,y,Nx+1,Ny+1,t(n));
    u(BD.Nord,n+1)=EDP.ubord(t(n+1),X(BD.Nord),Y(BD.Nord))';
    u(BD.Sud,n+1)=EDP.ubord(t(n+1),X(BD.Sud),Y(BD.Sud))';
    u(BD.Est,n+1)=EDP.ubord(t(n+1),X(BD.Est),Y(BD.Est))';
    u(BD.Ouest,n+1)=EDP.ubord(t(n+1),X(BD.Ouest),Y(BD.Ouest))';
  endfor
  
  
endfunction