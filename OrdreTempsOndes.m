close all

%ici on fixe le dx et on fait varier dt
N_pas=3;
pas=zeros(1,N_pas);
err=zeros(1,N_pas);

Nt=1000;
Nx=1000;


L=1; %longueur du domaine
T=2; %temps de la simulation

%Initialisation de la structure EDP
EDP.a=-L; EDP.b=L;
EDP.t0=0; EDP.T=T;
EDP.c=1; EDP.m=3;
EDP.uex=@(t,x) cos((2*EDP.m+1)*pi/(2*L)*x)*cos(EDP.c*(2*EDP.m+1)*pi/(2*L)*t);
EDP.u0=@(x) cos((2*EDP.m+1)*pi/(2*L)*x);
EDP.u1=@(x) 0*x;
EDP.ua=@(t) 0*t;
EDP.ub=@(t) 0*t;


for i=1:N_pas
  ht=EDP.T/Nt; hx=(EDP.b-EDP.a)/Nx;
  CFL=EDP.c*ht/hx;
  s=sprintf('Pour i=%d, CFL : si %f <= 1 => convergence',i,CFL);
  disp(s);
  
  [t,x,u]=EulerExplicite(EDP,Nt,Nx);
  
  pas(i)=ht;
  err(i)=max(max(abs(u-CalculF(EDP.uex,t,x)),[],2));
  Nt=2*Nt;
end

pente=(log(err(N_pas))-log(err(1)))/(log(pas(N_pas))-log(pas(1)));

figure(1);
loglog(pas,err,'r');
hold on;
loglog(pas,pas.^2,'ko-');
title("Représentation de l'erreur en fonction de dt");
legend(strcat('Erreur(dt), pente : ',num2str(pente)),"O(dt^2)","location", "southeast");
xlabel("dt");
ylabel('Erreur en norme L^\infty');
grid on;
