function TimeDomain_SS_Fraction(ModelID)
%% Prepare  Set of Reactions

global Kvec Uini DX
PrepareRxns

%% Initialize Variables To Store Parameters and Results of the Merged Model
EnsembleSize = 1000;
SSFrac = [0 0 0];
T_Rxn = 1e9;

for reps = 1:3

%% Analyze

SaveFileName = strcat('Time Domain Model%s Results',ModelID)
 
Xini = ones(length(X),1);
Uini = ones(length(U),1);

for Model = 1:EnsembleSize,

%Original
%     Kvec = (((ParamRange(:,1))-(ParamRange(:,2))).*rand(length(KVEC), 1)+(ParamRange(:,2)));
%     Kvec(ParamInfo(:,1)) = 10.^(-1 + 2*rand(size(ParamInfo,1),1));

%All uniform
    Kvec = (((ParamRange(:,1))-(ParamRange(:,2))).*rand(length(KVEC), 1)+(ParamRange(:,2)));
    Kvec(ParamInfo(:,1)) = 0.1 + 9.9*rand(size(ParamInfo,1),1);

%All log
%     Min = ParamRange(:,1);
%     Range = ParamRange(:,2) - ParamRange(:,1);
%     
%     Kvec = Min + 10.^(rand(length(KVEC) ,1).*log10(Range));
%     Kvec(ParamInfo(:,1)) = 10.^(-1 + 2*rand(size(ParamInfo,1),1));
    
%      Kvec(ParamInfo(7,1)) = Kvec(ParamInfo(8,1))/4;
%      Kvec(ParamInfo(8,1)) = Kvec(ParamInfo(8,1))*3/4;
    Kvec = Kvec';

    [Time,X_out]=ode15s('dxdt',[0 T_Rxn],Xini);
    Results{Model} = [Time X_out];
    if max(abs(X_out(end,:) - X_out(end - 1,:))) < 1e-3 && ...
        min(abs(X_out(end,:))) > 1e-6;
        SSFrac(reps) = SSFrac(reps) + 1
        Model
    end
end
end
disp(SSFrac/EnsembleSize)
%save(strcat(SaveFileName,'NoMetabRemove'),'Results');