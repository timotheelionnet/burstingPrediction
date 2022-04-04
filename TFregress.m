function [predScores,betaTrainAct,betaTrainInt]...
    = TFregress(M,Mnorm,tfList,charTFind,burstingData,A,I)

nTraining = numel(charTFind);
%predictIndex = setdiff(1:numel(tfList),charTFind);
predictIndex = 1:numel(tfList);
nPredict =numel(predictIndex);
%InteractionIndexStoreAct = betas_act_abs(1:300,2);
%InteractionIndexStoreInt = betas_int_abs(1:1000,2);
InteractionIndexStoreAct = find(sum(M)>=A);
InteractionIndexStoreInt = find(sum(M)>=I);
predScores = [];

[~,~,~,~,betaTrainAct,~,~,~] ...
    = plsregress(Mnorm(charTFind,InteractionIndexStoreAct),burstingData.activity,4);
predScores(:,1) = [ones(nPredict,1), Mnorm(predictIndex,InteractionIndexStoreAct)] * betaTrainAct;

[~,~,~,~,betaTrainInt,~,~,~] ...
    = plsregress(Mnorm(charTFind,InteractionIndexStoreInt),burstingData.intensity,6);
predScores(:,2) = [ones(nPredict,1), Mnorm(predictIndex,InteractionIndexStoreInt)] * betaTrainInt;

predTrainScoreact = [ones(nTraining,1), Mnorm(charTFind,InteractionIndexStoreAct)] * betaTrainAct;
predTrainScoreint = [ones(nTraining,1), Mnorm(charTFind,InteractionIndexStoreInt)] * betaTrainInt;

predScores(:,3) = 1:length(predScores);

end