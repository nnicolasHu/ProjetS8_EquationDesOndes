clear all;
close all;

L=1; %longueur du domaine
T=2; %temps de la simulation
  
Nt =100; Nx=35; Ny=35;

% 1. Initialisation de la structure EDP
EDP.a=0; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.k1=9; EDP.k2=3;
EDP.omega = EDP.c * sqrt(EDP.k1^2 + EDP.k2^2);
EDP.uex=@(t,x,y) cos(-EDP.omega*t + EDP.k1*x + EDP.k2*y);
EDP.u0=@(x,y) cos(EDP.k1*x + EDP.k2*y);
EDP.u1=@(x,y) EDP.omega * sin(EDP.k1*x + EDP.k2*y);
EDP.ubord=@(t,x,y) EDP.uex(t,x,y);
EDP.f =@(t,x,y) 0;




[t,x,y,u]=EulerExplicite2Dbis(EDP,Nt,Nx,Ny);
[yy xx] = meshgrid(y,x);

temps = 1.5;
indice = find((t-temps)==0);
for n=indice
  figure(1)
  surf(xx,yy,EDP.uex(t(n),xx,yy));
  title(strcat("Représentation de la solution exacte au temps t=", num2str(t(n)),", pour kx=", num2str(EDP.k1)));
  xlabel("x");
  ylabel("y");
  #surf(xx,yy,abs(EDP.uex(t(n),xx,yy)-Vec2dToMatrix(u(:,n),Nx+1,Ny+1)));
  figure(2)
  surf(xx,yy,Vec2dToMatrix(u(:,n),Nx+1,Ny+1));
  title(strcat("Représentation de la solution numérique au temps t=", num2str(t(n)),", pour kx=", num2str(EDP.k1)));
  xlabel("x");
  ylabel("y");
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