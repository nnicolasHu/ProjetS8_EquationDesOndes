close all

%ici on fixe le dt et on fait varier le dx
N_pas=5;
pas=zeros(1,N_pas);
err=zeros(1,N_pas);

Nt=2000;
Nx=100;


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
  
  pas(i)=hx;
  err(i)=max(max(abs(u-CalculF(EDP.uex,t,x))));
  Nx=2*Nx;
end

pente=(log(err(N_pas))-log(err(1)))/(log(pas(N_pas))-log(pas(1)));

figure(1);
loglog(pas,err,'r');
hold on;
loglog(pas,pas.^2,'ko-');
title("Représentation de l'erreur en fonction de h");
legend(strcat('Erreur(h), ordre2 : ',num2str(pente)),"O(h^2)","location", "southeast");
xlabel("h");
ylabel('Erreur en norme L^\infty');
grid on;




