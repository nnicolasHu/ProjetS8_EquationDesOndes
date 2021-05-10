clear all
close all

PLOT=1;
FREQ =10;

%constante du problème
m=3;
c=1;

%constante de discretisation
t0=0;
T=2;
L=1; %domaine en espace [-L,L]

Nx=100;
Nt=1000;

dt = (T-t0)/Nt;
dx = 2*L/Nx;

t = dt*[0:Nt];
x = dx*[0:Nx]-L;

uex=@(t,x) cos((2*m+1)*pi/(2*L)*x)*cos(c*(2*m+1)*pi/(2*L)*t);
U = CalculF(uex,t,x);
if PLOT
  figure(1);
  PlotSol(t,x,U,FREQ,'solution exacte',[x(1) x(end) -1 1],0.1);
end
