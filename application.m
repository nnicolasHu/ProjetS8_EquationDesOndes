clear all
close all


% 1.Initialisation de la structure EDP
L=6;
EDP.a=0;EDP.b=L;
EDP.t0=0;EDP.T=10;
nu=0.1;
% 2. Parametres de discretisation
Nx=500; 
Nt=5000;

disp('Calcul en cours...')
tic
%[t,x,u]=barre1(EDP,Nt,Nx,nu);
[t,x,u]=barre2(EDP,Nt,Nx);
tpsexc=toc
disp('Fin du calcul.')