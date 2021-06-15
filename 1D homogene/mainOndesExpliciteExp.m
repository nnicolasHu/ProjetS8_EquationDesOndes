clear all
close all

PLOT = 0;
FREQ = 10;
PLOT_EXACT = 1;

L=1; %longueur du domaine
T=2; %temps de la simulation

% 1. Initialisation de la structure EDP_reel
EDP_reel.a=-L; EDP_reel.b=L;
EDP_reel.t0=0; EDP_reel.T=T;
EDP_reel.c=1; EDP_reel.k=1;
EDP_reel.lambda=2*pi/EDP_reel.k;
EDP_reel.freq=EDP_reel.c/EDP_reel.lambda;
EDP_reel.omega=2*pi*EDP_reel.freq;
EDP_reel.uex=@(t,x) cos(EDP_reel.omega*t-EDP_reel.k*x);
EDP_reel.u0=@(x) cos(-EDP_reel.k*x);
EDP_reel.u1=@(x) EDP_reel.omega*cos(-EDP_reel.k*x);
EDP_reel.ua=@(t) cos(EDP_reel.omega*t-EDP_reel.k*EDP_reel.a);
EDP_reel.ub=@(t) cos(EDP_reel.omega*t-EDP_reel.k*EDP_reel.b);

% 2. Parametres de discretisation

Nx=1000;
Nt=1500;

% 3. Verification de la condition de C.F.L.
ht=EDP_reel.T/Nt; hx=(EDP_reel.b-EDP_reel.a)/Nx;
CFL=EDP_reel.c*ht/hx;
s=sprintf('CFL : si (%f <= 1) => convergence',CFL);
disp(s);

% 4. Resolution par le schema d'EULER explicite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Fonction non fournie : A IMPLEMENTER
disp('Calcul en cours...')
[t,x,u_relle]=EulerExplicite(EDP_reel,Nt,Nx);
disp('Fin du calcul.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5. Representations graphiques
Uex_relle=real(CalculF(EDP_reel.uex,t,x)); % Solution exacte
MIN=min(min(Uex_relle));
MAX=max(max(Uex_relle));
if PLOT
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

if PLOT
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
    
    %erreur
    figure();
    plot(x,abs(u_relle(:,i)-Uex_relle(:,i)));
    xlabel("x");
    title(strcat("Représentation de l'erreur en fonction de x en t=", num2str(t(i))));
  endfor
  
  
end