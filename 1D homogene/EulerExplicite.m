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
  u(:,2)=( spMatDiag(ones(1,Nx+1)) + (cte/2)*Lap1D(Nx+1))*EDP.u0(x)' + dt*EDP.u1(x)';
  
  u(1,:)=EDP.ua(t); %bord haud (x=-L)
  u(end,:)=EDP.ub(t); %bord bas (x=L)
  
  M = spMatDiag(2*ones(1,Nx+1))+ cte*Lap1D(Nx+1);
  for n=2:Nt
    u(:,n+1) = M*u(:,n) - u(:,n-1);
    u(1,n+1) = EDP.ua(t(n+1));
    u(end,n+1) = EDP.ub(t(n+1));
  endfor
  
  
endfunction