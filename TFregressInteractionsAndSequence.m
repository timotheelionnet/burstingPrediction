if exist('ActData','var') == 0
    PLS_regression_interactionANDsequence; close all
end
%%
A=270;
I=40;
M3act = M3norm(:,sum(M3)>=A);
M3int = M3norm(:,sum(M3)>=I);
ActData = [M3act,parametersSubAct];
IntData = [M3int,parametersSubInt];

ncomp = [2,3];
trainingIndex = 1:height(burstingData2);
predictIndex = 1:numel(tfList2);
nPredict =numel(predictIndex);
tic
%%
predScores = [];

% train model for active fraction, return weights for the interactors
[XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
            = plsregress(ActData(trainingIndex,:),burstingData2.activity(trainingIndex),ncomp(1));
% predict active fraction for all TFs
predScores(:,1) = [ones(nPredict,1), ActData(predictIndex,:)] * betaTrainAct;

% train model for intensity
[XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
            = plsregress(IntData(trainingIndex,:),burstingData2.intensity(trainingIndex),ncomp(2));
% predict intensity for all TFs
predScores(:,2) = [ones(nPredict,1), IntData(predictIndex,:)] * betaTrainInt;

predScores(:,3) = 1:length(predScores);

% from biogridScript
plotPredictions(predScores,burstingData2,tfList2,...
charTFind,exampleTFInd,characterizedInd)
        