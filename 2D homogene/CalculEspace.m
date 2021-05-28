function F=CalculEspace(f,x,y)
%commentaire à refaire
% Calcul de la fonction f(t,x) 
% INPUT :
%   f :  fonction du type @(t,x)
%   t :  discretisation en temps, vecteur de dimension nt
%   x : discretisation en espace,  vecteur de dimension nx
% OUTPUT :
%   F : matrice definie par F(i,n)=f(t(n),x(i)), de dimension (nx,nt)
ny=length(y);
nx=length(x);
F=zeros(ny,nx);
for n=1:nx
    F(:,n)=f(x(n),y)';
end