clear all;
close all;

L=1; %longueur du domaine
T=2; %temps de la simulation

##N_pas=2;
##pas=zeros(1,N_pas);
##err=zeros(1,N_pas);
Nt =21; Nx=10; Ny=10;

% 1. Initialisation de la structure EDP
EDP.a=0; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1;
EDP.uex=@(t,x,y) sin(pi*x).*sin(pi*y).*cos(t);
EDP.u0=@(x,y) sin(pi*x).*sin(pi*y);
EDP.u1=@(x,y) zeros(1,length(x));
EDP.ubord=@(t,x,y) 0;
EDP.f =@(t,x,y) (EDP.c^2 *2* pi^2 -1)*sin(pi*x).*sin(pi*y).*cos(t);



[t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny);
[yy xx] = meshgrid(y,x);

for n=15
  figure(1)
  surf(xx,yy,EDP.uex(t(n),xx,yy));
  #surf(xx,yy,abs(EDP.uex(t(n),xx,yy)-Vec2dToMatrix(u(:,n),Nx+1,Ny+1)));
  figure(2)
  surf(xx,yy,Vec2dToMatrix(u(:,n),Nx+1,Ny+1))
endfor
##for i=1:N_pas
##  hy=(EDP.t0-EDP.T)/Ny;
##  [t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny);
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
##  disp(i)
##  pas(i)=hy;
##  abs(abs(u)-abs(Uex))
##  err(i)=max(max(abs(abs(u)-abs(Uex))));
##  Ny=2*Ny;
##end

##err
##pente=(log(err(N_pas))-log(err(1)))/(log(pas(N_pas))-log(pas(1)));
##
##figure(1);
##loglog(pas,err,'r');
##hold on;
##loglog(pas,pas.^2,'ko-');
##title("Représentation de l'erreur en fonction de dx");
##legend(strcat('Erreur(dx), pente : ',num2str(pente)),"O(h^2)","location", "southeast");
##xlabel("dx");
##grid on;
##ylabel('Erreur en norme L^\infty');