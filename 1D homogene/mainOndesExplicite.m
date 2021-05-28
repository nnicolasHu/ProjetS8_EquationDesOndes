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
EDP.c=1; EDP.m=3;
EDP.uex=@(t,x) cos((2*EDP.m+1)*pi*x/(2*L))*cos(EDP.c*(2*EDP.m+1)*pi*t/(2*L));
EDP.u0=@(x) cos((2*EDP.m+1)*pi*x/(2*L));
EDP.u1=@(x) 0*x;
EDP.ua=@(t) 0*t;
EDP.ub=@(t) 0*t;

% 2. Parametres de discretisation

Nx=2000;
Nt=2000;

% 3. Verification de la condition de C.F.L.
ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
CFL=EDP.c*ht/hx;
s=sprintf('CFL : si (%f <= 1) => convergence',CFL);
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
MIN=min(min(Uex));
MAX=max(max(Uex));
if PLOT
    figure(1)
    PlotSol2(t,x,u,'freq',FREQ,'title','Sol. Appr.','axis',[x(1) x(end) MIN MAX],'pause',0.1)
    % PlotSol(t,x,u,FREQ,'Sol. Appr.',[x(1) x(end) MIN MAX],0.1) % Old version
end

Err=abs(u-Uex)/MAX; % relative error
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

fprintf('Erreur relative  (max en temps et espace) : %e\n',max(Ninf)/max(max(abs(Uex))));

% 6. Représentation graphique pour certains temps
if PLOT_EXACT
  for i=[2 ceil(Nt/2) Nt+1]
    figure();
    plot(x,u(:,i),";solution numérique;");
    hold on;
    plot (x,Uex(:,i),";solution exacte;");
    legend();
    xlabel("x");
    title(strcat("Représentation de la solution en fonction de x en t=", num2str(t(i))));
    hold off;
    
    %erreur
    figure();
    plot(x,abs(u(:,i)-Uex(:,i)));
    xlabel("x");
    title(strcat("Représentation de l'erreur en fonction de x en t=", num2str(t(i))));
  endfor
  
  
end

