<<<<<<< HEAD
function F=CalculF(f,t,x)
% Calcul de la fonction f(t,x) 
% INPUT :
%   f :  fonction du type @(t,x)
%   t :  discretisation en temps, vecteur de dimension nt
%   x : discretisation en espace,  vecteur de dimension nx
% OUTPUT :
%   F : matrice definie par F(i,n)=f(t(n),x(i)), de dimension (nx,nt)
nt=length(t);
nx=length(x);
F=zeros(nx,nt);
for n=1:nt
    F(:,n)=f(t(n),x)';
=======
function F=CalculF(f,t,x)
% Calcul de la fonction f(t,x) 
% INPUT :
%   f :  fonction du type @(t,x)
%   t :  discretisation en temps, vecteur de dimension nt
%   x : discretisation en espace,  vecteur de dimension nx
% OUTPUT :
%   F : matrice definie par F(i,n)=f(t(n),x(i)), de dimension (nx,nt)
nt=length(t);
nx=length(x);
F=zeros(nx,nt);
for n=1:nt
    F(:,n)=f(t(n),x)';
>>>>>>> de207de9f8c99e922567b870b4ca49fe8cd3d992
end