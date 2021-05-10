clear all
close all

PLOT = 1;
FREQ = 10;

L=1 %longueur du domaine
T=2 %temps de la simulation

% 1. Initialisation de la structure EDP
EDP.a=-L; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.m=3;
EDP.uex=@(t,x) cos((2*EDP.m+1)*pi/(2*L)*x)*cos(EDP.c*(2*EDP.m+1)*pi/(2*L)*t);
EDP.u0=@(x) cos((2*EDP.m+1)*pi/(2*L)*x);
EDP.u1=@(x) 0*x;


% 2. Parametres de discretisation

Nx=500;
Nt=2100;

% 3. Verification de la condition de C.F.L.

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

