clear all
close all

PLOT_EXACT = 1;

L=1; %longueur du domaine
T=2; %temps de la simulation

% 1. Initialisation de la structure EDP
EDP.a=0; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.k1=[3,6,9]; EDP.k2=[3,6,9];


% 6. Représentation graphique pour certains temps
Nx=10;Ny=10;
for indice=1:length(EDP.k1)
  k1 = EDP.k1(indice);
  k2 = 3;
  EDP.omega = EDP.c * sqrt(k1^2 + k2^2);
  EDP.uex=@(t,x,y) cos(-EDP.omega*t + k1*x + k2*y);
  EDP.u0=@(x,y) cos(k1*x + k2*y);
  EDP.u1=@(x,y) EDP.omega * sin(k1*x + k2*y);
  EDP.ubord=@(t,x,y) EDP.uex(t,x,y);
  EDP.f =@(t,x,y) 0;
##  EDP.uex=@(t,x,y) cos(k1*x).*cos(k2*y).*cos(t);
##  EDP.u0=@(x,y) cos(k1*x).*cos(k2*y);
##  EDP.u1=@(x,y) 0;
##  EDP.ubord=@(t,x,y) EDP.uex(t,x,y);
##  EDP.f =@(t,x,y) (-1+k1^2+k2^2)*cos(k1*x).*cos(k2*y).*cos(t);
  Nt=100;
##  EDP.omega=EDP.c*sqrt(k1^2+k2^2);
##  EDP.uex=@(t,x,y) cos(EDP.omega*t-k1*x-k2*y);
##  EDP.u0=@(x,y) cos(-k1*x-k2*y);
##  EDP.u1=@(x,y) -EDP.omega*sin(-k1*x-k2*y);
##  EDP.ubord=@(t,x,y) EDP.uex(t,x,y);
##  EDP.f =@(t,x,y) 0;
##  Nt=100;
  err=1;
  ht=EDP.T/Nt;
  while (err > 10^(-1.5))
    Nx+=5;
    Ny=Nx;
    h=(EDP.b-EDP.a)/Nx;
    CFL=EDP.c*ht/h;
    s=sprintf('CFL : si (%f < %f) => convergence',CFL,sqrt(2)/2);
    disp(s);
    
    % 4. Resolution par le schema d'EULER explicite
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp('Calcul en cours...')
    [t,x,y,u]=EulerExplicite2Dbis(EDP,Nt,Nx,Ny);
    %disp('Fin du calcul.')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 5. Representations graphiques
    N=(Nx+1)*(Ny+1);
  
    Uex=zeros(N,Nt+1);
    [I,J]=bijRecF(1:N,Nx+1);
    X=x(I+1);
    Y=y(J+1);
    for j=1:Nt+1
      Uex(:,j)=(EDP.uex(t(j),X,Y));
    endfor
    MAX=max(max(Uex));
   
    Err=abs(u-Uex)/MAX; % relative error
    err=max(NormInf(Err))/max(max(abs(Uex)))
    
  endwhile
  
  fprintf('Erreur relative  (max en temps et espace) : %e\n',err);
  
##  figure();
##  plot(x,u(:,end),";solution numérique;");
##  hold on;
##  plot (x,Uex(:,end),";solution exacte;");
##  legend();
##  xlabel("x");
##  title(strcat("Représentation de la solution en fonction de x en t=", num2str(t(end))));
##  hold off;
##  
  s=sprintf('Pour k1=%d, il a fallu %d points de discrétisation en espace',k1,Nx);
  disp(s);
  err=1;
endfor
##      temps = 1.5;
##      indice = find((t-temps)==0);