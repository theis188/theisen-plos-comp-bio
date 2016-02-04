clear
figure
hold on
pinkish_rgb = [1,1-(1-0.811764705882353)/2,1-(1-0.862745098039216)/2;];
slate_rgb = 1-(1- [0.317647058823529,0.396078431372549,0.447058823529412;])/1.3

fill([1.0 1.0 1.45 1.45], [-4.5 3 3 -4.5], pinkish_rgb, 'EdgeColor', 'none')
plot(linspace(0,1),sin(exp(1.8*linspace(0,1))))
plot(1,sin(exp(1.8)),'bo','MarkerFaceColor','b')
%Xs
% plot(0.9,0*sin(exp(1.8*0.9))-4.5,'x','color','black','linewidth',2, ...
% 'MarkerSize', 10)
% plot(1.1,-4.5,'x','color','black','linewidth',2, ...
% 'MarkerSize', 10)
xlim([-0.2 1.4])
ylim([-12 3])

set(gca,'box', 'off')
set(gcf,'color','white' ,'position', [230 250 500 550])
axis off
%Axes
annotation('textarrow',[0.2 0.2],[0.465 0.85],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 2, ...
      'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)

annotation('textarrow',[0.15 0.95],[0.515 0.515],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 2, ...
      'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)
%Style
% annotation('textarrow',[0.76 0.76],[0.68 0.54],'string',' ', ...
%       'HeadStyle','none','LineStyle', ':', 'linewidth', 2)
%   
% annotation('textarrow',[0.76 0.76],[0.535 0.53],'string',' ', ...
%       'HeadStyle','vback2','LineStyle', '-', 'linewidth', 0.5)
%   
% annotation('textarrow',[0.662 0.662],[0.68 0.54],'string',' ', ...
%       'HeadStyle','none','LineStyle', '-', 'linewidth', 2)
%   
% annotation('textarrow',[0.662 0.662],[0.535 0.53],'string',' ', ...
%       'HeadStyle','vback2','LineStyle', '-', 'linewidth', 0.5)



text(1.06,2.5,'Unstable','fontname', 'Calibri', 'fontsize', 14)
text(1.1,1.9,'Region','fontname', 'Calibri', 'fontsize', 14)

text(0.33,2.5,'Stable','fontname', 'Calibri', 'fontsize', 14)
text(0.32,1.9,'Region','fontname', 'Calibri', 'fontsize', 14)
text(0.32,1.9,'Region','fontname', 'Calibri', 'fontsize', 14)

%line([0.9-0.131; 1.031-0.131], [0.0-2.5; 1.6*2/3-2.5], ...
%    'color','black','linewidth',2.5)
%line([0.9+0.07; 1.031+0.07], [0.0-6; 1.6*2/3-6], ...
%    'color','black','linewidth',2.5,'linestyle',':')