function Ninf=NormInf(U)
% Calcul de la norme L^\infty de U en espace
% INPUT :
%   U: matrice de dimension (nx,nt)
% OUTPUT :
%   Ninf :  vecteur de dimension nt
[nx,nt]=size(U);
Ninf=zeros(1,nt);
for n=1:nt
    Ninf(n)=norm(U(:,n),inf);
end