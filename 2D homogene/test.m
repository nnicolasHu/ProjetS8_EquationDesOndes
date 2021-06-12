clear all;
close all;

L=1; %longueur du domaine
T=2; %temps de la simulation

Nt =10; Nx=40; Ny=30;

% 1. Initialisation de la structure EDP
EDP.a=0; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1;
EDP.uex=@(t,x,y) sin(pi*x).*sin(pi*y).*cos(t);
EDP.u0=@(x,y) sin(pi*x).*sin(pi*y);
EDP.u1=@(x,y) 0;
EDP.ubord=@(t,x,y) 0;
EDP.f =@(t,x,y) (2*pi^2 -1) sin(pi*x).*sin(pi*y).*cos(t)


[t,x,y,u]=EulerExplicite2D(EDP,Nt,Nx,Ny);