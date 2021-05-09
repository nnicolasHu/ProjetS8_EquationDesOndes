function [t,x,u]=OndeEulerExplicite(EDP,Nt,Nx)
  t0=EDP.t0;
  T=EDP.T;
  a=EDP.a;
  b=EDP.b;
  dt=T/Nt;
  dx=(b-a)/Nx;
  t=t0:dt:T+t0;
  x=a:dx:b;
  c=EDP.c;
  alpha=((c*dt)/dx)^2;
  u=zeros(Nx+1,Nt+1);
  u1=zeros(Nx+1,1);
  u(:,1)=EDP.u0h(x);
  u1(:,1)=EDP.u1h(x);
  u(:,2)=(eye(Nx+1)-(dt*dt/2)*Lap1D(Nx+1))*u(:,1)+dt*u1(:,1);
  for j=2:Nt
    for i=2:Nx
      u(i,j+1)=2*(1-alpha^2)*u(i,j)+(alpha^2)*(u(i-1,j)+u(i+1,j))-u(i,j-1);
    endfor
  endfor
endfunction
