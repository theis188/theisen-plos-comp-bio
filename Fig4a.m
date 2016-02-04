clear
load('FluxFigData.mat')

DBColor = [0.3765 0.4863 0.5569];
DRColor = [0.7333 0.2471 0.2471];

U = ones(11,1)
Uf1 = U
Uf1([7 8]) = 10
Uf2 = U
Uf2([7 8]) = 0.1

for i = 1:5
    for PertStep = 1:StepsUp + 1
        OverallFluxes = FLUXES(GoodResults(:,PertStep,i), ...
        EnsembleKvec(:,i), 1, U + (Uf1 - U)*(PertStep - 1/200));
        FpkFlux(i,PertStep) = OverallFluxes(11,1)/2;
    end
    RightBound(i,:) = 200 + find(isnan(FpkFlux(i,:)),1)
end

for i = 1:5
    for PertStep = 1:StepsUp + 1
        OverallFluxes = FLUXES(GoodResults(:,201 + PertStep,i), ...
        EnsembleKvec(:,i), 1, U + (Uf2 - U)*(PertStep - 1/200));
        FpkFluxDown(i,PertStep) = OverallFluxes(11,1)/2;
    end
    LeftBound(i,:) = 201 - find(isnan(FpkFluxDown(i,:)),1)
end

XAXIS = [log10(linspace(0.1,1,201)) log10(linspace(1,10,201))];
FpkFlux = [fliplr(FpkFluxDown) FpkFlux];

figure
hold on

 

for i = 1
    plot(XAXIS(LeftBound(i):RightBound(i)),...
        FpkFlux(i,LeftBound(i):RightBound(i)), ...
        'LineWidth', 2+0.5*(i==1), 'color', DBColor)
    plot(XAXIS(RightBound(i)), ...
        FpkFlux(i,RightBound(i)), 'color', DBColor, 'markersize', 8 + (i==1)*6,...
         'marker', 'pentagram', 'markerfacecolor', DBColor)
end
%     plot([XAXIS(find(isnan(FpkFlux(1,:)),1)-1) log10(1.6)],...
%         [FpkFlux(1,find(isnan(FpkFlux(1,:)),1) - 1)... 
%         FpkFlux(1,find(isnan(FpkFlux(1,:)),1) - 1)], ...
%         'r:', 'LineWidth', 2)
    %plot(log10(1.6), FpkFlux(1,find(isnan(FpkFlux(1,:)),1) - 1), ...
    %    'ro', 'markersize', 8, 'markerfacecolor', 'red')
    
annotation('textarrow',[0.742 0.742],[0.445 0.15], 'linewidth', 2, ...
      'HeadStyle','none','LineStyle', '--', 'Color', DRColor);  

  annotation('textarrow',[0.742 0.742],[0.445 0.135], 'linewidth', 2, ...
      'HeadStyle','vback2', 'headlength', 12 ,'headwidth', 12 ,'LineStyle', 'none', 'Color', DRColor);
   
%plot(XAXIS(find(isnan(FpkFlux(9,:)),1)-1), ...
%    FpkFlux(9,find(isnan(FpkFlux(9,:)),1) - 1), 'r', 'markersize', 10, 'marker', 'pentagram', 'markerfacecolor', 'r')
set(gca,'LineWidth',1)
set(gca,'xtick',[-log10(2) 0 log10(2) 5 10], 'xticklabel',{'0.5' '1' '2' '5' '10'}, 'fontsize', 16, 'fontname', 'Cambria')
set(gca,'ytick',[0 0.5 1], 'yticklabel',{'0' '0.5' '1'})
xlabel('Fold Change Phosphoketolase', 'fontsize', 16, 'fontname', 'Cambria')
ylabel('Acetate Flux')
ylim([0 1.1])
xlim([-log10(2) log10(2)])
set(gcf,'color','white','position', [230 250 500 380])