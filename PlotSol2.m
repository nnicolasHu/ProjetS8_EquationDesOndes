function PlotSol2(t,x,u,varargin)
% PlotSol2(t,x,u)
%   representation au cours des temps t d'une donnée u (temps/espace) 
%   t est la discretisation en temps, x la discrétisation en espace et
%   u(n,i) la valeur au temps t(n) et au point x(i).
% Cette fonction accepte des arguments optionnels sous la forme de paires
% de clef/valeur
%  'freq': pour specifier la frequence d'affichage (ie 1:freq:length(t))
%            PlotSol2(t,x,u,'freq',3)
%          Defaut: 1
%  'pause': pour spécifier la durée d'une pause entre chaque représentation
%            PlotSol2(t,x,u,'pause',0.1)
%          Defaut: 0 (pas de pause)
%  'axis': pour spécifier l'axe de représentation 
%            PlotSol2(t,x,u,'axis',[0 10 -1 1])
%          Defaut: calculer automatiquement
%  'title': pour spécifier le titre affiché
%            PlotSol2(t,x,u,'title','essai')
%  Bien evidemment, on peut spécifier plusieurs clef/valeurs :
%    PlotSol2(t,x,u,'title','essai','freq',5)
  p = inputParser;
  p.addParamValue('freq',1,@(x) isscalar(x) && x>=1);
  p.addParamValue('pause',0,@isscalar);
  p.addParamValue('axis',[],@(x) isempty(x) || (size(x,1)==1 && size(x,2)==4));
  p.addParamValue('title','',@ischar);
  p.parse(varargin{:});
  R=p.Results;
  
  nt=length(t);
  format=''; % create string format for title
  if length(R.title)>0
    format=[R.title,' - temps t=%f'];
  end
  R
  for n=1:R.freq:nt
      plot(x,u(:,n))
      if ~isempty(R.axis)
        axis(R.axis)
      end
      title(sprintf(format,t(n)));
      xlabel('x')
      drawnow
      if (R.pause >0), pause(R.pause); end
  end
  
end
