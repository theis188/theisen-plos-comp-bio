%declare global
clear
global Kvec Uini DX
load('MCC_EquationModel')





%% Analyze
NoEnzymes = length(EnzName);

Xini = ones(length(X),1);
Uini = ones(length(U),1);
Results = {};

T_Rxn = 100;
%Vm_min = 0.5;
%Vm_max = 1.5;
%Keq = 1;
%Km_min = 0.1;
%Km_max = 1;
%log10_min = -4;
%log10_max = 1;
%steps = 100;
%EnzChange = 'Xpk';
%EnzInd=find(ismember(EnzName,EnzChange));

Fpk_Locs = [];
Xpk_Locs = [];
for i=1:size(EnzName,1);
    switch EnzName{i,1};
        case 'Fpk';
            Fpk_Locs = [Fpk_Locs i];
        case 'Xpk';
            Xpk_Locs = [Xpk_Locs i];
    end
end

FinalMetabs = Net.MetabName
%FinalMetabs(ToRemove, :) = []
figure
hold on

DB = [0 0.0118 0.3569]
DR = [0.5176 0 0]
SemiDB = mean([[1 1 1]; DB; DB])
SemiDR = mean([[1 1 1]; DR; DR])

% for U_step = [1 2];
for U_step = [1 1.1 1.7 1.8 2.0 10];


Uini(Fpk_Locs,1) = U_step;
Uini(Xpk_Locs,1) = U_step;
    
Xini(find(ismember(FinalMetabs,'MeO')),1) = 200;
% Xini(find(ismember(FinalMetabs,'R5P')),1) = 2;
% %Xini(find(ismember(FinalMetabs,'NAD')),1) = 2

%Kvec = [30.5772175435402;4.31126460340185;3.24674027892026;1.99865345322414;370.127342461046;4.01527658881462;9.70969128775782;4.61063647294729;4.76319456670636;18.6826603894600;8.31152526561847;0.778540310238755;0.334112468388380;118.833483379548;4.38661689278484;4.78209888679055;5.45509293687614;6.45020030547977;6.69073144061450;128.310200019946;6.48638318443191;7.63049493514619;5.18590730695563;3.06907331442723;8.64102371799638;115.008455321247;1.62952102422410;1.79556798180984;5.19896764879040;3.46828428686665;5.48162928905627;5.33023350310018;9.66046700620036;15.8558061731322;9.57053744875481;107.712738394267;6.52827038739583;4.71282249095874;0.118134114281791;8.25760637968777;9.86720051483171;2.27290220498415;5.19563444088776;2;]';
Kvec = [75.9913161758704;4.85241539051879;6.01912442580630;0.746458422783380;1110.38555885090;1.47057391941565;9.70177656854063;5.63314299059581;3.06297941627980;31.6393304243756;5.40850101987372;3.93878551870894;2.61094745491610;39.7656771260256;6.96947519162648;4.41547923966235;1.53286590290923;4.16732982718948;4.96841215254077;148.797444767071;3.84581492143534;8.94725126507607;5.03062738986026;9.81025001214123;5.61816710644470;37.3633128220707;9.59797338514058;0.642611037700055;8.21423443906270;6.33790935422621;7.96813308471008;3.35793237274586;5.71586474549172;15.0022555832132;9.00150372214214;20.6299596313430;7.84869333827134;6.53894682127778;4.47320332608857;19.4911537768977;6.20853462512486;6.18355330560464;6.23142858612071;2]';
Kvec([1 44]) = 0;

    [Time,X_out]=ode15s('dxdt',[0 T_Rxn],Xini);
    %[uf, conc1, Bif1] = SimpleODESolver(StepsUp, Xini, U, Uf1, Kvec,DVDX, DVDU,S,JACOBIAN );

Results{end + 1} = [Time X_out];
DxDt = [];
for t = 1:size(Time,1) - 1;
    DxDt(end + 1,:) = (X_out(t+1,:) - X_out(t,:))/(Time(t+1) - Time(t));
end


color = DB;
plot(Results{1,size(Results,2)}(:,1),Results{1,size(Results,2)}(:,11), ...
   'Color', [1 1 1] - ([1 1 1] - color) * (7 - size(Results,2)) / 6, ...
   'LineWidth', 3)



% color = DR;
% plot(Results{1,size(Results,2)}(1:end-1,1),DxDt(:,10), ...
%     'Color', [1 1 1] - ([1 1 1] - color) * (7 - size(Results,2)) / 6, ...
%     'LineWidth', 3)

end

% legend({'1x PK', ...
%     '1.7x PK', ...
%     '1.8x PK', ...
%     '2x PK', ...
%     '2.4x PK', ...
%     '10x PK'}, ...
%     'FontSize', 12, 'Position',[0.2 0.7 0.2 0.1]);
% legend('boxoff');





% figure 
% hold on
% plot(Results{1,1}(:,1),Results{1,1}(:,9),'LineStyle', '-', 'color',SemiDB,'LineWidth',3)
% plot(Results{1,1}(:,1),Results{1,1}(:,10),'LineStyle', '--', 'color',SemiDB,'LineWidth',3)
% plot(Results{1,2}(:,1),Results{1,2}(:,9),'LineStyle', '-', 'color',SemiDR,'LineWidth',3)
% plot(Results{1,2}(:,1),Results{1,2}(:,10),'LineStyle', '--', 'color',SemiDR,'LineWidth',3)
% 
set(gca,'xtick',[0 50 100], 'xticklabel',{'0' '50' '100'}, 'fontsize', 16, 'fontname', 'Cambria')
set(gca,'ytick',[0 60 120], 'yticklabel',{'0' '60' '120'})
xlabel('Time', 'fontsize', 16, 'fontname', 'Cambria')
ylabel('Concentration')
ylim([0 120])
xlim([0 100])
set(gcf,'color','white','position', [230 250 500 380])

% legend({'X5P, 1x PK', ...
%     'R5P, 1x PK', ...
%     'X5P, 2x PK', ...
%     'R5P, 2x PK'}, ...
%     'FontSize', 12, 'Position',[0.2 0.7 0.2 0.1]);
% legend('boxoff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(Results{1,1}(:,1),Results{1,1}(:,7),'LineStyle', '-', 'color',SemiDB,'LineWidth',3)
% plot(Results{1,1}(:,1),Results{1,1}(:,5),'LineStyle', '--', 'color',SemiDB,'LineWidth',3)
% plot(Results{1,2}(:,1),Results{1,2}(:,7),'LineStyle', '-', 'color',SemiDR,'LineWidth',3)
% plot(Results{1,2}(:,1),Results{1,2}(:,5),'LineStyle', '--', 'color',SemiDR,'LineWidth',3)
% 
% set(gca,'xtick',[0 50 100], 'xticklabel',{'0' '50' '100'}, 'fontsize', 16, 'fontname', 'Cambria')
% set(gca,'ytick',[0 4 8], 'yticklabel',{'0' '4' '8'})
% xlabel('Time', 'fontsize', 16, 'fontname', 'Cambria')
% ylabel('Concentration')
% ylim([0 8])
% xlim([0 100])
% set(gcf,'color','white','position', [230 250 500 380])
% 
% legend({'G3P, 1x PK', ...
%     'F6P, 1x PK', ...
%     'G3P, 2x PK', ...
%     'F6P, 2x PK'}, ...
%     'FontSize', 12, 'Position',[0.3 0.45 0.2 0.1]);
% legend('boxoff');

 

% for n = 1:6
%     X_Out(n) = Results{n}(end,11)
% end
% 
% for i = 1:6
%     h = bar(i, X_Out(i))
%     BColor = [1 1 1] - ([1 1 1] - (color)) * (7 - i) / 6
%     set(h, 'FaceColor', BColor)
% end
% 
% set(gca,'xtick',[1 2 3 4 5 6],'xticklabel',{'1x' '1.7x' '1.8x' '2x' '2.4x' '10x'}, ...
%     'fontsize', 16, 'fontname', 'Cambria')
% set(gca,'ytick',[0 60 120], 'yticklabel',{'0' '60' '120'})
% xlabel('PK Fold Change', 'fontsize', 16, 'fontname', 'Cambria')
% ylabel('Final Acetate Titer')
% ylim([0 120])
% xlim([0 7])
% set(gcf,'color','white','position', [230 250 500 380])

