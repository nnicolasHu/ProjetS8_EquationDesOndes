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
end