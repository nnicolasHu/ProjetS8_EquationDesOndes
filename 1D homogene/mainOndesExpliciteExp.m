clear all
close all

PLOT = 1;
FREQ = 15;
PLOT_EXACT = 0;

L=1; %longueur du domaine
T=2; %temps de la simulation

% 1. Initialisation de la structure EDP
EDP.a=-L; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.k=2;
EDP.lambda=2*pi/EDP.k;
EDP.freq=EDP.c/EDP.lambda;
EDP.omega=2*pi*EDP.freq;
EDP.uex=@(t,x) cos(EDP.omega*t-EDP.k*x);
EDP.u0=@(x) cos(-EDP.k*x);
EDP.u1=@(x) -EDP.omega*sin(-EDP.k*x);
EDP.ua=@(t) cos(EDP.omega*t-EDP.k*EDP.a);
EDP.ub=@(t) cos(EDP.omega*t-EDP.k*EDP.b);

% 2. Parametres de discretisation

Nx=100;
Nt=1000;

% 3. Verification de la condition de C.F.L.
ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
CFL=EDP.c*ht/hx;
s=sprintf('CFL : si (%f <= 1) => convergence',CFL);
disp(s);

% 4. Resolution par le schema d'EULER explicite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Fonction non fournie : A IMPLEMENTER
disp('Calcul en cours...')
[t,x,u_relle]=EulerExplicite(EDP,Nt,Nx);
disp('Fin du calcul.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5. Representations graphiques
Uex_relle=CalculF(EDP.uex,t,x); % Solution exacte
MIN=min(min(Uex_relle));
MAX=max(max(Uex_relle));
if PLOT && 0
    figure(1)
    PlotSol2(t,x,u_relle,'freq',FREQ,'title','Sol. Appr.','axis',[x(1) x(end) MIN MAX],'pause',0.1)
    % PlotSol(t,x,u,FREQ,'Sol. Appr.',[x(1) x(end) MIN MAX],0.1) % Old version
end

Err=abs(u_relle-Uex_relle)/MAX; % relative error
MAX=max(max(Err));
if MAX >1
    AXIS=[];
else 
    AXIS=[x(1) x(end) 0 MAX];
end

if PLOT && 0
    figure(2)
    % PlotSol(t,x,Err,FREQ,'Erreur',AXIS,0.1)
    PlotSol2(t,x,Err,'freq',FREQ,'title','Erreur','axis',AXIS,'pause',0.1)
end

Ninf=NormInf(Err);
if PLOT
  figure(3)
  plot(t,Ninf);
  xlabel('t')
  title('Erreur en norme L^\infty en espace')
end

fprintf('Erreur relative  (max en temps et espace) : %e\n',max(Ninf)/max(max(abs(Uex_relle))));

% 6. Représentation graphique pour certains temps
if PLOT_EXACT
  for i=[2 ceil(Nt/2) Nt+1]
    figure();
    plot(x,u_relle(:,i),";solution numérique;");
    hold on;
    plot (x,Uex_relle(:,i),";solution exacte;");
    legend();
    xlabel("x");
    title(strcat("Représentation de la solution en fonction de x en t=", num2str(t(i))));
    hold off;
    
  endfor
  
  
end