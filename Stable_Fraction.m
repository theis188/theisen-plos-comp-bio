function Stable_Fraction(ModelID)
%% Prepare  Set of Reactions
PrepareRxns_ConservedGroups

%% Initialize Variables To Store Parameters and Results of the Merged Model
EnsembleSize = 1e3;
% EnsembleKvec = NaN(length(KVEC), EnsembleSize);
        
%% Analyze
NoEnzymes = length(EnzName);

Xini = ones(length(X),1);
Uini = ones(length(U),1);

TotStab = 0

for Model = 1:EnsembleSize,
    Stable = 0;
    Kvec = (((ParamRange(:,1))-(ParamRange(:,2))).*rand(length(KVEC), 1)+(ParamRange(:,2)));
    Kvec(ParamInfo(:,1)) = K1S(Xini, Kvec, 1, rVnet,ones(length(rVnet),1));
    Stable = max(real(eig(JACOBIAN(0, Xini, Kvec, 1,ones(length(rVnet),1)))))<-1e-10;
    TotStab = TotStab + Stable;
end

fprintf('Stable Fraction is: %f \n',TotStab/EnsembleSize)