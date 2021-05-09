function PlotSol(t,x,u,Freq,titre,AXIS,tpause)
% INPUT :
%     Freq : Frequence de representation graphique en temps (voir boucle for)
%   tpause : Temps de pause entre 2 courbes en temps si tpause>0,
%               pas de pause sinon
%     AXIS : utilisé pour la commande axis
% COMMENTAIRES A FINIR
nt=length(t);

for n=1:Freq:nt
    plot(x,u(:,n))
    if ~isempty(AXIS)
       axis(AXIS)
    end
    s=sprintf('%s - temps t=%f',titre,t(n));
    title(s)
    xlabel('x')
    drawnow
    if (tpause >0)
        pause(tpause)
    end
end