close all


Nt=1000;
Nx=1000;

delta=1/10;
L=1;
T=2;
c0=1;
dt=T/Nt;
dx=2*L/Nx;
x=-L:dx:L;
t=0:dt:T;
w=@(x) heaviside(x+delta)*heaviside(delta-x);
u=@(t,x) 1/2*(w(x-c0*t)+w(x+c0*t))-1/2*(w(x+2*L-c0*t)+ w(x-2*L+c0*t))+1/2*(w(x+4*L-c0*t)+w(x-4*L+c0*t));

Uex=CalculF(u,t,x); % Solution exacte
MIN=min(min(Uex));
MAX=max(max(Uex));
if PLOT

    PlotSol2(t,x,Uex,'freq',FREQ,'title','Sol.Exacte Propagation Créneau.','axis',[x(1) x(end) MIN MAX],'pause',0.1)
end