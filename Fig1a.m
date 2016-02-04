clear
figure
hold on
plot([1.0 2.75], [0.8 -0.6], 'linewidth',2)
plot([2],[0],'bo','MarkerFaceColor','b')
ylim([-2.5 2.5])
xlim([-0.5 4])

set(gca,'box', 'off')
set(gcf,'color','white' ,'position', [230 250 400 530])
axis off
%Axes
annotation('textarrow',[0.32 0.32],[0.5 0.72],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 2, ...
      'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)

annotation('textarrow',[0.32 0.32],[0.5 0.42],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 2, ...
      'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)
  
annotation('textarrow',[0.32 0.76],[0.515 0.515],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 2, ...
      'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)
%Small Arrows  
annotation('textarrow',[0.45 0.51],[0.615 0.57],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 0.5, ...
      'HeadLength', 6, 'HeadWidth', 5, 'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)

annotation('textarrow',[0.515 0.55],[0.565 0.54],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 0.5, ...
      'HeadLength', 6, 'HeadWidth', 5, 'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)

annotation('textarrow',[0.69 0.625],[0.435 0.48],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 0.5, ...
      'HeadLength', 6, 'HeadWidth', 5, 'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)

annotation('textarrow',[0.62 0.585],[0.485 0.51],'string',' ', ...
      'HeadStyle','vback2','LineStyle', '-', 'linewidth', 0.5, ...
      'HeadLength', 6, 'HeadWidth', 5, 'TextRotation',90, ...
      'fontname', 'Cambria', 'fontsize', 16)
  


%Labels
text(2.8,-0.2,'X','FontName','Cambria', 'fontsize', 16)

text(0.25,0.02,'0','FontName','Cambria', 'fontsize', 16)
text(0.05,1.1,'dX','FontName','Cambria', 'fontsize', 16)
text(0.07,0.76,'dt','FontName','Cambria', 'fontsize', 16)

line([0.05; 0.4], [0.91; 0.91], 'color', 'black','linewidth', 2)

