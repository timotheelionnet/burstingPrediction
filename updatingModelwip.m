function [predScores,betaTrainAct,betaTrainInt] = updatingModelwip(M,Mnorm,tfList,charTFind,burstingData,A,I)

nComp = 10;
trainingIndex = charTFind;
nTraining = numel(trainingIndex);
validIndex = setdiff(1:numel(tfList),trainingIndex);
nValid =numel(validIndex);
%InteractionIndexStoreAct = betas_act_abs(1:300,2);
%InteractionIndexStoreInt = betas_int_abs(1:1000,2);
InteractionIndexStoreAct = find(sum(M)>=A);
InteractionIndexStoreInt = find(sum(M)>=I);
predScores = [];

[XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,MSEtrainact,statstrainact] ...
= plsregress(Mnorm(trainingIndex,InteractionIndexStoreAct),burstingData.activity,nComp);
predScores(:,1) = [ones(nValid,1), Mnorm(validIndex,InteractionIndexStoreAct)] * betaTrainAct;

[XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,MSEtrainint,statstrainint] ...
= plsregress(Mnorm(trainingIndex,InteractionIndexStoreInt),burstingData.intensity,nComp);
predScores(:,2) = [ones(nValid,1), Mnorm(validIndex,InteractionIndexStoreInt)] * betaTrainInt;

predTrainScoreact = [ones(nTraining,1), Mnorm(trainingIndex,InteractionIndexStoreAct)] * betaTrainAct;
predTrainTraincorrAct = corr(burstingData.activity,predTrainScoreact);

predTrainScoreint = [ones(nTraining,1), Mnorm(trainingIndex,InteractionIndexStoreInt)] * betaTrainInt;
predTrainTraincorrInt = corr(burstingData.intensity,predTrainScoreint);

predScores(:,3) = 1:length(predScores);

end