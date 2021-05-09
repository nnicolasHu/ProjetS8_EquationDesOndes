clear all
close all

PLOT = 1;
FREQ = 10;

% 1. Initialisation de la structure EDP
EDP.a=0;EDP.b=2*pi;
EDP.t0=0;EDP.T=2;
EDP.c=2;
k=5;
EDP.uex=@(t,x) cos(k*t)*cos(x);
EDP.uprim=@(t,x) -k*sin(k*t)*cos(x);
EDP.u0h=@(x) EDP.uex(EDP.t0,x);
EDP.u1h=@(x) EDP.uprim(EDP.t0,x);

% 2. Parametres de discretisation

Nx=100;
%Nt=2000; % => divergence
Nt=2100;% => convergence 

%Nx=150;
%Nt=4500; % => divergence
%Nt=5000; % => convergence 

%  Nx=200;
% Nt = 8000; % => divergence
% Nt = 8400; % => convergence

% 3. Verification de la condition de C.F.L.
ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
CFL=EDP.c*ht/(hx^2);
s=sprintf('CFL : si %f < 0.5 => convergence',CFL);
disp(s);

% 4. Resolution par le schema d'EULER explicite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Fonction non fournie : A IMPLEMENTER
disp('Calcul en cours...')
%[t,x,u]=EulerExpliciteMix(EDP,Nt,Nx,1,0,1,0);
%[t,x,u]=EulerExplicitebis(EDP,Nt,Nx);
[t,x,u]=OndeEulerExplicite(EDP,Nt,Nx);
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

