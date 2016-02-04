function Main(ModelID)
%% Prepare  Set of Reactions
PrepareRxns_ConservedGroups

%% Initialize Variables To Store Parameters and Results of the Merged Model
EnsembleSize = 1000;
% EnsembleKvec = NaN(length(KVEC), EnsembleSize);
        
%% Analyze
NoEnzymes = length(EnzName);
StepsUp = 200;
StepsDown = 200;
PertUp = 10;
PertDown = .1;

ModelResults = repmat([{NaN(size(S,1),(1+StepsUp)+(1+StepsDown),NoEnzymes)} ...
    {NaN(NoEnzymes,1)} {NaN(NoEnzymes,1)}], EnsembleSize, 1);

SaveFileName = sprintf('Model%s Results 1000',ModelID)

Xini = ones(length(X),1);
Uini = ones(length(U),1);



for Model = 1:EnsembleSize,
    Model
    Stable = 0;
    while ~Stable,
        Kvec = (((ParamRange(:,1))-(ParamRange(:,2))).*rand(length(KVEC), 1)+(ParamRange(:,2)));
        Kvec(ParamInfo(:,1)) = K1S(Xini, Kvec, 1, rVnet,ones(length(rVnet),1));
        Stable = max(real(eig(JACOBIAN(0, Xini, Kvec, 1,ones(length(rVnet),1)))))<-1e-10;
        %real(eig(JACOBIAN(0, Xini, Kvec, 1,ones(length(rVnet),1))))
        %pause
    end
    Model
        
    EnsembleKvec(:,Model) = Kvec;
    
    EnzymeConc = NaN(size(S,1),(1+StepsUp)+(1+StepsDown),NoEnzymes);
    EnzymeBifur = NaN(NoEnzymes,1);
    ElapsedTime = NaN(NoEnzymes,1);
        
    for Enzyme=1:NoEnzymes,

        TimeIn = clock;
        U = Uini;
        Uf1 = U;
        Uf1(Enzyme) = PertUp;
        Uf2 = U;
        Uf2(Enzyme) = PertDown;

        [uf, conc1, Bif1] = SimpleODESolver(StepsUp, Xini, U, Uf1, Kvec,DVDX, DVDU,S,JACOBIAN );
        [uf2, conc2, Bif2] = SimpleODESolver(StepsDown, Xini, U, Uf2, Kvec,DVDX, DVDU,S,JACOBIAN );
        
        ElapsedTime(Enzyme) = etime(clock,TimeIn);
        EnzymeConc(:,:,Enzyme) = [conc1',conc2'];
        EnzymeBifur(Enzyme) = Bif1+2*Bif2;
        
    end
    
    ModelResults(Model,:) = [{EnzymeConc} {EnzymeBifur} {ElapsedTime}];

end

save(SaveFileName);