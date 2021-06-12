clear all;
close all;

L=1; %longueur du domaine
T=2; %temps de la simulation

Nt =10; Nx=40; Ny=30;

% 1. Initialisation de la structure EDP
EDP.a=0; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.k=2;
EDP.lambda=2*pi/EDP.k;
EDP.f=EDP.c/EDP.lambda;
EDP.omega=2*pi*EDP.f;
EDP.uex=@(t,x,y) sin(pi*x).*sin(pi*y).*cos(t);
EDP.u0=@(x,y) sin(pi*x).*sin(pi*y);
EDP.u1=@(x) i*EDP.omega*exp(-i*EDP.k*x);
EDP.ua=@(t) exp(i*(EDP.omega*t-EDP.k*EDP.a));
EDP.ub=@(t) exp(i*(EDP.omega*t-EDP.k*EDP.b));


[t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny);