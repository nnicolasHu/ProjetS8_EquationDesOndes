clear all
close all

PLOT = 1;
FREQ = 10;

% 1. Initialisation de la structure EDP
m=3;
l=1;
EDP.a=-l;EDP.b=l;
EDP.t0=0;EDP.T=2;
EDP.c=1;
EDP.uex=@(t,x) 0.5*(cos((2*m+1)*pi*(x+EDP.c*t)/(2*l))+cos((2*m+1)*pi*(x-EDP.c*t)/(2*l)));
EDP.u0=@(x) cos(((2*m+1)*pi*x)/(2*l));
EDP.ua=@(t) 0*t;
EDP.ub=@(t) 0*t;
EDP.u1=@(x) 0*x;

dx=zeros(2,1);
Err=zeros(2,1);
% 2. Parametres de discretisation


Nx=500; % => divergence
Nt=2100;% => convergence 
%Nx=[100 200 400];
%Nx=150;
%Nt=4500; % => divergence
%Nt=5000; % => convergence 

%  Nx=200;
% Nt = 8000; % => divergence
% Nt = 8400; % => convergence
%k=1;
%for i=Nx
  
  % 3. Verification de la condition de C.F.L.
  ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/i;
  CFL=EDP.c*ht/hx;
  s=sprintf('CFL : si %f < 1 => convergence',CFL);
  disp(s);

  % 4. Resolution par le schema d'EULER explicite
  
  disp('Calcul en cours...')
  [t,x,u,dx(k),dt]=OndeEulerExplicite(EDP,Nt,Nx);
  disp('Fin du calcul.')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Uex=CalculF(EDP.uex,t,x);
%  Err(k)=max(max(abs(u-Uex)));
%  k=k+1;
%end



%dx
%Err
%loglog(dx,Err);
%hold on
%loglog(dx,dx.^2,'ks-');
% 5. Representations graphiques
Uex=CalculF(EDP.uex,t,x); % Solution exacte
MIN=min(min(Uex));
MAX=max(max(Uex));
if PLOT
    figure(1)
    PlotSol2(t,x,u,'freq',FREQ,'title','Sol. Appr.','axis',[x(1) x(end) MIN MAX],'pause',0.1)
   
    PlotSol2(t,x,Uex,'freq',FREQ,'title','Sol. Ex.','axis',[x(1) x(end) MIN MAX],'pause',0.1)
    %PlotSol(t,x,u,FREQ,'Sol. Appr.',[x(1) x(end) MIN MAX],0.1) % Old version
%end




Err=abs(u-Uex)/MAX; % relative error
MAX=max(max(Err));
if MAX >1
    AXIS=[];
else 
    AXIS=[x(1) x(end) 0 MAX];
end

%if PLOT
%    figure(2)
    % PlotSol(t,x,Err,FREQ,'Erreur',AXIS,0.1)
%    PlotSol2(t,x,Err,'freq',FREQ,'title','Erreur','axis',AXIS,'pause',0.1)
%end

Ninf=NormInf(Err);
%if PLOT
 % figure(3)
 % plot(t,Ninf);
 % xlabel('t')
 % title('Erreur en norme L^\infty en espace')
%end

fprintf('Erreur relative  (max en temps et espace) : %e\n',max(Ninf)/max(max(abs(Uex))));