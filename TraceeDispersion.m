close all

alpha = [0.2:0.2:1];
G = linspace(0,0.5,50);

%calcul de q
q = 1./(pi*alpha'*G) .* asin(alpha'*sin(pi*G));

labels = {};
figure();
hold on;
for i=1:size(q,1)
  plot(G,q(i,:));
  labels = {labels{:}, ["alpha : ", num2str(alpha(i))]};
endfor
title("Courbes de dispersion pour différents alpha");
xlabel("G");
ylim([min(min(q))-0.02,1.02]);
legend(labels,"location","southwest");
