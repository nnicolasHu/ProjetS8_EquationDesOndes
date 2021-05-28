close all

%ici on fixe le dx et on fait varier dt
N_pas=2;
pas=zeros(1,N_pas);
err=zeros(1,N_pas);

Nt=2000;
Nx=1800;


L=1; %longueur du domaine
T=2; %temps de la simulation

%Initialisation de la structure EDP
EDP.a=-L; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.m=3;
EDP.uex=@(t,x) cos((2*EDP.m+1)*pi*x/(2*L))*cos(EDP.c*(2*EDP.m+1)*pi*t/(2*L));
EDP.u0=@(x) cos((2*EDP.m+1)*pi*x/(2*L));
EDP.u1=@(x) 0*x;
EDP.ua=@(t) 0*t;
EDP.ub=@(t) 0*t;


for i=1:N_pas
  ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
  CFL=EDP.c*ht/hx;
  s=sprintf('Pour i=%d, CFL : si %f < 1 => convergence',i,CFL);
  disp(s);
  
  [t,x,u]=EulerExplicite(EDP,Nt,Nx);
  
  pas(i)=ht;
  err(i)=max(max(abs(u-CalculF(EDP.uex,t,x))));
  Nt=2*Nt;
end

pente=(log(err(N_pas))-log(err(1)))/(log(pas(N_pas))-log(pas(1)));

figure(1);
loglog(pas,err,'r');
hold on;
loglog(pas,1./pas.^2,'ko-');
title("Représentation de l'erreur en fonction de dt");
legend(strcat('Erreur(dt), pente : ',num2str(pente)),"O(1/dt^2)","location", "southwest");
xlabel("dt");
ylabel('Erreur en norme L^\infty');
grid on;
