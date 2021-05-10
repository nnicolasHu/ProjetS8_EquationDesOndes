function [t,x,u]=EulerExplicite(EDP,Nt,Nx)
  dt=(EDP.T-EDP.t0)/Nt;
  dx=(EDP.b-EDP.a)/Nx;
  
  t=EDP.t0 + dt*[0:Nt];
  x=EDP.a + dx*[0:Nx];
  
  cte=(EDP.c*dt/dx)**2;
  
  %initialisation
  u=zeros(Nx+1,Nt+1); %colonne=espace et ligne=temps
  u(:,1)=EDP.u0(x)'; %bord gauche (t=t0)
  %Calcul de u en t=dt (u(:,2))
  Ah=Lap1D(Nx+1);
  Ah(1,Nx)=1;
  Ah(Nx+1,2);
  
  u(:,2)=( spMatDiag(ones(1,Nx+1)) + cte/2*Ah)*EDP.u0(x)' + dt*EDP.u1(x)';
  
  for n=2:Nt
    u(:,n+1) = (spMatDiag(2*ones(1,Nx+1))+ cte*Ah)*u(:,n) - u(:,n-1);  
  endfor
  
  
endfunction