function [predScores,betaTrainAct,betaTrainInt,InteractorIndAct,InteractorIndInt] ...
            = TFregress(M,Mnorm,tfList,charTFind,burstingData,...
                        minInteractionsA,minInteractionsI,nCompA,nCompI)
% this functions trains a regression model to predict the bursting kinetics (active fraction and intensity)
% from interaction data, and predicts the kinetics of all TFs in the
% interaction matrix.

%%%%% INPUT
% M: interaction matrix, with row indices corresponding to gene symbols in
        % tfList.
% Mnorm: same as M, but z-score transformed.       
% tfList: list of gene symbols corresponding to the rows of M.
% charTFind: list of row indices corresponding to the TFs that have been
        % experimentally characterized.
% burstingData: MATLAB table with columns TFnames, activity, intensity (CRISPRburst results)  
% A: minimum number of interactions required for an interactor (column of
        % M) to make in order to be retained in the regression to predict active
        % fraction
% I: same as A, but for intensity regression

%%%%% OUTPUT
% predScores: array with size <number of TFs> x 3 where:
        % col 1 is the predicted active fraction for each TF, 
        % col 2 is the predicted intensity for each TF
        % col 3 is the TF index corresponding to the gene symbol list tfList
        
% betaTrainAct: predictive weight for each of the interactors retained in
        % the active fraction model, relative to the index list 
        % in InteractorIndAct
% betaTrainInt: predictive weight for each of the interactors retained in
        % the intensity model, relative to the index list in
        % InteractorIndInt
% InteractorIndAct: the indices of the columns of M retained in the model
        % for active fraction
% InteractorIndInt: the indices of the columns of M retained in the model
        % for intensity
%%
predictIndex = 1:numel(tfList);
nPredict =numel(predictIndex);

% collect indices of interactors that pass the threshold for each model.
InteractorIndAct = find(sum(M)>=minInteractionsA);
InteractorIndInt = find(sum(M)>=minInteractionsI);

predScores = [];

% train model for active fraction, return weights for the interactors
[~,~,~,~,betaTrainAct,~,~,~] ...
    = plsregress(Mnorm(charTFind,InteractorIndAct),burstingData.activity,nCompA);

% predict active fraction for all TFs
predScores(:,1) = [ones(nPredict,1), Mnorm(predictIndex,InteractorIndAct)] * betaTrainAct;

% train model for intensity
[~,~,~,~,betaTrainInt,~,~,~] ...
    = plsregress(Mnorm(charTFind,InteractorIndInt),burstingData.intensity,nCompI);

% predict intensity for all TFs
predScores(:,2) = [ones(nPredict,1), Mnorm(predictIndex,InteractorIndInt)] * betaTrainInt;

predScores(:,3) = 1:length(predScores);

end