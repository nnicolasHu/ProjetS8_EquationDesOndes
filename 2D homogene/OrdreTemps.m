clear all;
close all;

L=1; %longueur du domaine
T=4; %temps de la simulation

N_pas=3;
pas=zeros(1,N_pas);
err=zeros(1,N_pas);
tmp = 16;
Nt =96; 
Nx=tmp; Ny=tmp;

% 1. Initialisation de la structure EDP
EDP.a=0; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1;
EDP.uex=@(t,x,y) sin(pi*x).*sin(pi*y).*cos(t^2);
EDP.u0=@(x,y) sin(pi*x).*sin(pi*y);
EDP.u1=@(x,y) 0;
EDP.ubord=@(t,x,y) EDP.uex(t,x,y);
EDP.f =@(t,x,y) sin(pi*x).*sin(pi*y).*(-2.*sin(t.^2)-4.*t.^2.*cos(t.^2)) + EDP.c^2 *2* pi^2*sin(pi*x).*sin(pi*y).*cos(t.^2) ;

h=(EDP.b-EDP.a)/Nx
dt=(EDP.T-EDP.t0)/Nt

[t,x,y,u]=EulerExplicite2Dbis(EDP,Nt,Nx,Ny);
[yy xx] = meshgrid(y,x);

for n=30
  figure(1)
  surf(xx,yy,EDP.uex(t(n),xx,yy));
  %surf(xx,yy,abs(EDP.uex(t(n),xx,yy)-Vec2dToMatrix(u(:,n),Nx+1,Ny+1)));
  figure(2)
  surf(xx,yy,Vec2dToMatrix(u(:,n),Nx+1,Ny+1))
endfor

for i=1:N_pas
  h=(EDP.b-EDP.a)/Ny;
  ht=EDP.T/Nt;
  CFL=EDP.c*ht/h;
  s=sprintf('Pour i=%d, CFL : si %f < %f => convergence',i,CFL,sqrt(2)/2);
  disp(s);
  [t,x,y,u]=EulerExplicite2Dbis(EDP,Nt,Nx,Ny);
  
  N=(Nx+1)*(Ny+1);
  
  Uex=zeros(N,Nt+1);
  [I,J]=bijRecF(1:N,Nx+1);
  X=x(I+1);
  Y=y(J+1);
  for j=1:Nt+1
    Uex(:,j)=(EDP.uex(t(j),X,Y));
  endfor
  pas(i)=ht;
  err(i)=max(max(abs(u-Uex)));
  Nt=Nt+17;
end

err
pente=(log(err(N_pas))-log(err(1)))/(log(pas(N_pas))-log(pas(1)));

figure();
loglog(pas,err,'r');
hold on;
loglog(pas,pas.^2,'ko-');
xlim([0.025,0.05]);
title("Representation de l'erreur en temps");
legend(strcat('Erreur(dx), pente : ',num2str(pente)),"O(dt^2)","location", "southeast");
xlabel("dt");
grid on;
ylabel('Erreur en norme L^\infty');


##clear all;
##close all;
##
##L=1; %longueur du domaine
##T=2; %temps de la simulation
##
##N_pas=2;
##pas=zeros(1,N_pas);
##err=zeros(1,N_pas);
##Nt =30; Nx=10; Ny=10;
##
##% 1. Initialisation de la structure EDP
##EDP.a=0; EDP.b=L;
##EDP.t0=0; EDP.T=T;
##EDP.c=1;
##EDP.uex=@(t,x,y) y.*exp(x)+t; % sin(pi*x).*sin(pi*y).*cos(t);
##EDP.u0=@(x,y) y.*exp(x); %sin(pi*x).*sin(pi*y);
##EDP.u1=@(x,y) 1;
##EDP.ubord=@(t,x,y) EDP.uex(t,x,y);
##EDP.f =@(t,x,y) -EDP.c^2*(y.*exp(x)); %(EDP.c^2 *2* pi^2 -1)*sin(pi*x).*sin(pi*y).*cos(t);
##
##h=(EDP.b-EDP.a)/Nx;
##dt=(EDP.T-EDP.t0)/Nt;
##
##[t,x,y,u]=EulerExplicite2Dbis(EDP,Nt,Nx,Ny);
##[yy xx] = meshgrid(y,x);
##
##%for n=15
##%  figure(1)
##%  surf(xx,yy,EDP.uex(t(n),xx,yy));
##  #surf(xx,yy,abs(EDP.uex(t(n),xx,yy)-Vec2dToMatrix(u(:,n),Nx+1,Ny+1)));
##%  figure(2)
##%  surf(xx,yy,Vec2dToMatrix(u(:,n),Nx+1,Ny+1))
##%endfor
##
##for i=1:N_pas
##  h=(EDP.b-EDP.a)/Ny;
##  ht=EDP.T/Nt;
##  CFL=EDP.c*ht/h;
##  s=sprintf('Pour i=%d, CFL : si %f < %f => convergence',i,CFL,sqrt(2)/2);
##  disp(s);
##  [t,x,y,u]=EulerExplicite2Dbis(EDP,Nt,Nx,Ny);
##  
##  N=(Nx+1)*(Ny+1);
##  
##  Uex=zeros(N,Nt+1);
##  [I,J]=bijRecF(1:N,Nx+1);
##  X=x(I+1);
##  Y=y(J+1);
##  for j=1:Nt+1
##    Uex(:,j)=(EDP.uex(t(j),X,Y));
##  endfor
##  pas(i)=ht;
##  err(i)=max(max(abs(u-Uex)));
##  Nt=Nt+20;
##end
##
##err
##pente=(log(err(N_pas))-log(err(1)))/(log(pas(N_pas))-log(pas(1)));
##
##figure(3);
##loglog(pas,err,'r');
##hold on;
##loglog(pas,pas.^2,'ko-');
##title("Repr?sentation de l'erreur en temps");
##legend(strcat('Erreur(dx), pente : ',num2str(pente)),"O(h^2)","location", "southeast");
##xlabel("dt");
##grid on;
##ylabel('Erreur en norme L^\infty');