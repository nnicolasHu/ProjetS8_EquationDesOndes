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
Nx=50;
% 6. Représentation graphique pour certains temps
if PLOT_EXACT
  for m=EDP.m
     Nt=1000;
     err=1;
     ht=EDP.T/Nt;
     EDP.uex=@(t,x) cos((2*m+1)*pi*x/(2*L))*cos(EDP.c*(2*m+1)*pi*t/(2*L));
     EDP.u0=@(x) cos((2*m+1)*pi*x/(2*L));
     EDP.u1=@(x) 0*x;
     EDP.ua=@(t) 0*t;
     EDP.ub=@(t) 0*t;
     while (err > 10^(-2))
      Nx+=50;
      hx=(EDP.b-EDP.a)/Nx;
      CFL=EDP.c*ht/hx;
      s=sprintf('CFL : si (%f < 1) => convergence',CFL);
      disp(s);
      s=sprintf('CFL : si (%f < 1) => convergence',CFL);
      disp(s);
      

      
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
      plot(x,u(:,end),";solution numérique;");
      hold on;
      plot (x,Uex(:,end),";solution exacte;");
      legend();
      xlabel("x");
      title(strcat("Représentation de la solution en fonction de x en t=", num2str(t(end))));
      hold off;
      
      Nx=Nx+15;
      ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
      CFL=EDP.c*ht/hx;
      %erreur
      %figure();
      %plot(x,abs(u(:,i)-Uex(:,i)));
      %xlabel("x");
      %title(strcat("Représentation de l'erreur en fonction de x en t=", num2str(t(i))));
    endwhile
   s=sprintf('Pour k=%d, il a fallu %d points de discrétisation en espace',m,Nx);
   disp(s);
   err=1;;
  endfor
  
  
end