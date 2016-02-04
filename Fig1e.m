clear
figure
hold on
pinkish_rgb = [1,1-(1-0.811764705882353)/2,1-(1-0.862745098039216)/2;];
slate_rgb = 1-(1- [0.317647058823529,0.396078431372549,0.447058823529412;])/1.3

X = linspace(0,1)
Y = sin(exp(1.8*X))
for i = 1:size(X,2)-1
DYDX(i) = (Y(i+1) - Y(i)) / (X(i+1) - X(i))
end
plot(X(1:end-1),DYDX)
%Xs
xlim([-0.2 1.4])
ylim([-25 15])

set(gca,'box', 'off')
set(gcf,'color','white' ,'position', [230 250 500 550])
axis off
%Axes
annotation('textarrow',[0.2 0.2],[0.465 0.85],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 2, ...
      'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)

annotation('textarrow',[0.15 0.95],[0.61 0.61],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 2, ...
      'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)
%Style  



%line([0.9-0.131; 1.031-0.131], [0.0-2.5; 1.6*2/3-2.5], ...
%    'color','black','linewidth',2.5)
%line([0.9+0.07; 1.031+0.07], [0.0-6; 1.6*2/3-6], ...
%    'color','black','linewidth',2.5,'linestyle',':')