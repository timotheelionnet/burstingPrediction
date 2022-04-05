function [predScores,betaTrainAct,betaTrainInt] = updatingModelwip(M,Mnorm,tfList,charTFind,burstingData,A,I)

nComp = 10;
nTraining = numel(charTFind);
predictIndex = setdiff(1:numel(tfList),charTFind);
nPredict =numel(predictIndex);
%InteractionIndexStoreAct = betas_act_abs(1:300,2);
%InteractionIndexStoreInt = betas_int_abs(1:1000,2);
InteractionIndexStoreAct = find(sum(M)>=A);
InteractionIndexStoreInt = find(sum(M)>=I);
predScores = [];

[XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,MSEtrainact,statstrainact] ...
= plsregress(Mnorm(charTFind,InteractionIndexStoreAct),burstingData.activity,nComp);
predScores(:,1) = [ones(nPredict,1), Mnorm(predictIndex,InteractionIndexStoreAct)] * betaTrainAct;

[XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,MSEtrainint,statstrainint] ...
= plsregress(Mnorm(charTFind,InteractionIndexStoreInt),burstingData.intensity,nComp);
predScores(:,2) = [ones(nPredict,1), Mnorm(predictIndex,InteractionIndexStoreInt)] * betaTrainInt;

predTrainScoreact = [ones(nTraining,1), Mnorm(charTFind,InteractionIndexStoreAct)] * betaTrainAct;
predTrainTraincorrAct = corr(burstingData.activity,predTrainScoreact);

predTrainScoreint = [ones(nTraining,1), Mnorm(charTFind,InteractionIndexStoreInt)] * betaTrainInt;
predTrainTraincorrInt = corr(burstingData.intensity,predTrainScoreint);

predScores(:,3) = 1:length(predScores);

end