function [t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny)
  dt=(EDP.T-EDP.t0)/Nt;
  dx=(EDP.b-EDP.a)/Nx;
  dy=(EDP.b-EDP.a)/Ny;
  
  t=EDP.t0 + dt*[0:Nt];
  x=EDP.a + dx*[0:Nx];
  y=EDP.a + dy*[0:Ny];
  
  %Bord
  BD = BordDomaine(Nx+1,Ny+1);
  idxBD=unique([BD.Sud,BD.Nord,BD.Est,BD.Ouest]);
  idxBDc=setdiff(1:(Ny+1)*(Nx+1),idxBD);
  Ac=sparse(idxBD,idxBD,ones(1,length(idxBD)),(Ny+1)*(Nx+1),(Ny+1)*(Nx+1));

    %initialisation en espace
  [I,J]=bijRecF(1:(Ny+1)*(Nx+1),Nx+1);
  X=x(I+1);
  Y=y(J+1);
  u(BD.Nord,:)=EDP.ubord(t,X(BD.Nord),Y(BD.Nord));
  u(BD.Sud,:)=EDP.ubord(t,X(BD.Sud),Y(BD.Sud));
  u(BD.Est,:)=EDP.ubord(t,X(BD.Est),Y(BD.Est));
  u(BD.Ouest,:)=EDP.ubord(t,X(BD.Ouest),Y(BD.Ouest));
  

  %initialisation en temps
  u=zeros((Ny+1)*(Nx+1),Nt+1); %le carré en espace est "compresse" en un vecteur colone
  u(:,1) = EvalFun2DVec(EDP.u0,x,y,Nx+1,Ny+1);
  Lap2D = Lap2DAssembling(Nx+1,Ny+1,dx,dy)+Ac;
  u(:,2) = u(:,1) + dt*EvalFun2DVec(EDP.u1,x,y,Nx+1,Ny+1)' + ( dt^2 /2)*(EvalSource2DVec(EDP.f,x,y,Nx+1,Ny+1,t(1)) + EDP.c^2*Lap2D*u(:,1));
 
  for n=2:Nt
    u(:,n+1) = ((EDP.c*dt)^2)*Lap2D*u(:,n) + 2*u(:,n) - u(:,n-1)+ (dt^2)*EvalSource2DVec(EDP.f,x,y,Nx+1,Ny+1,t(n));
  endfor
  
  
endfunction