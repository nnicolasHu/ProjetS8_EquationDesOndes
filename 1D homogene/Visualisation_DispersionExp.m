clear all
close all

PLOT_EXACT = 1;

L=1; %longueur du domaine
T=2; %temps de la simulation

% 1. Initialisation de la structure EDP
EDP.a=-L; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.k=[10,20,40];


% 6. Représentation graphique pour certains temps
Nx=50;
for indice=1:length(EDP.k)
  k = EDP.k(indice);
  EDP.lambda=2*pi/k;
  EDP.freq=EDP.c/EDP.lambda;
  EDP.omega=2*pi*EDP.freq;
  EDP.uex=@(t,x) cos(EDP.omega*t-k*x);
  EDP.u0=@(x) cos(-k*x);
  EDP.u1=@(x) -EDP.omega*sin(-k*x);
  EDP.ua=@(t) cos(EDP.omega*t-k*EDP.a);
  EDP.ub=@(t) cos(EDP.omega*t-k*EDP.b);
  
  Nt=1000;
  err=1;
  ht=EDP.T/Nt;
  while (err > 10^(-2))
    Nx+=50;
    hx=(EDP.b-EDP.a)/Nx;
    CFL=EDP.c*ht/hx;
    s=sprintf('CFL : si (%f < 1) => convergence',CFL);
    disp(s);
    
    % 4. Resolution par le schema d'EULER explicite
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp('Calcul en cours...')
    [t,x,u]=EulerExplicite(EDP,Nt,Nx);
    %disp('Fin du calcul.')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 5. Representations graphiques
    Uex=CalculF(EDP.uex,t,x); % Solution exacte
    MAX=max(max(Uex));
   
    Err=abs(u-Uex)/MAX; % relative error
    err=max(NormInf(Err))/max(max(abs(Uex)));
    
  endwhile
  
  fprintf('Erreur relative  (max en temps et espace) : %e\n',err);
  
  figure();
  plot(x,u(:,end),";solution numérique;");
  hold on;
  plot (x,Uex(:,end),";solution exacte;");
  legend();
  xlabel("x");
  title(strcat("Représentation de la solution en fonction de x en t=", num2str(t(end))));
  hold off;
  
  s=sprintf('Pour k=%d, il a fallu %d points de discrétisation en espace',k,Nx);
  disp(s);
  err=1;
endfor
##      temps = 1.5;
##      indice = find((t-temps)==0);