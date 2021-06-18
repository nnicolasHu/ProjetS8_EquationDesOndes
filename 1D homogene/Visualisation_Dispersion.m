clear all
close all

PLOT = 0;
FREQ = 10;
PLOT_EXACT = 1;

L=1; %longueur du domaine
T=2; %temps de la simulation

% 1. Initialisation de la structure EDP
EDP.a=-L; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.m=[3,6,9];

% 6. Représentation graphique pour certains temps
if PLOT_EXACT
  for m=EDP.m
     Nx=500;
     Nt=100000;
     err=1;
     ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
     CFL=EDP.c*ht/hx;
     while (err>10^(-2)) & (CFL<1)
      s=sprintf('CFL : si (%f < 1) => convergence',CFL);
      disp(s);
      EDP.uex=@(t,x) cos((2*m+1)*pi*x/(2*L))*cos(EDP.c*(2*m+1)*pi*t/(2*L));
      EDP.u0=@(x) cos((2*m+1)*pi*x/(2*L));
      EDP.u1=@(x) 0*x;
      EDP.ua=@(t) 0*t;
      EDP.ub=@(t) 0*t;

      
      % 4. Resolution par le schema d'EULER explicite
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %     Fonction non fournie : A IMPLEMENTER
      disp('Calcul en cours...')
      [t,x,u]=EulerExplicite(EDP,Nt,Nx);
      disp('Fin du calcul.')
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % 5. Representations graphiques
      Uex=CalculF(EDP.uex,t,x); % Solution exacte
      MAX=max(max(Uex));
   
      Err=abs(u-Uex)/MAX; % relative error
      err=max(NormInf(Err))/max(max(abs(Uex)));

      fprintf('Erreur relative  (max en temps et espace) : %e\n',err);
      

      figure();
      plot(x,u(:,2),";solution numérique;");
      hold on;
      plot (x,Uex(:,2),";solution exacte;");
      legend();
      xlabel("x");
      title(strcat("Représentation de la solution en fonction de x en t=", num2str(t(2))));
      hold off;
      
      Nx=Nx*2;
      ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
      CFL=EDP.c*ht/hx;
      %erreur
      %figure();
      %plot(x,abs(u(:,i)-Uex(:,i)));
      %xlabel("x");
      %title(strcat("Représentation de l'erreur en fonction de x en t=", num2str(t(i))));
    endwhile
   err=1;
   s=sprintf('Pour m=%d, il a fallu %d points de discrétisation en espace',m,Nx/2);
   disp(s);
  endfor
  
  
end