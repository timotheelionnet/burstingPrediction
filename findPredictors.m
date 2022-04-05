function [topPredictorsActivity, topPredictorsIntensity] ...
    = findPredictors(betaTrainAct, betaTrainInt, M, intList, minInteractionsA, minInteractionsI)

%%%%% INPUT

% betaTrainAct: predictive weight for each of the interactors retained in
        % the active fraction model. 
        % !!! Indices are relative to the columns of M, but once 
        % they've been filtered out for interactors that do not have enough
        % interactions. So make sure that the values of minInteractionsA 
        % and minInteractionsI fed to this
        % function match those used in the function that generates the
        % beta weights (TFregress), otherwise interactor names will be off.
        
% betaTrainInt: same as betaTrainAct, but for intensity. Same warning applies.

% M: interaction matrix, with rows corresponding to TFs and columns to
        % interactors. Column indices map to gene symbols in intList.
        
% intList: list of gene symbols corresponding to the columns of M.

% minInteractionsA: minimum number of interactions required for an interactor (column of
        % M) to make in order to be retained in the regression to predict active
        % fraction. !!!!! Make sure the value entered in this function
        % matches the value of the same parameter fed to the function that
        % generated the beta weights (TFregress).
        
% minInteractionsI: same as A, but for intensity regression. Same warning
        % applies.

%%%%% OUTPUT
% topPredictorsActivity: list of gene symbols corresponding to 100 interactors 
        % most highly predictive of the active fraction, based on the
        % weights beta.

% topPredictorsIntensity: list of gene symbols corresponding to 100 interactors 
        % most highly predictive of the intensity, based on the
        % weights beta.

%%  
% number of interactors to include in each list:
nInteractorsToOutput = 100;        
        
%% active fraction
InteractorIndAct = find(sum(M)>=minInteractionsA); % retain only interactors that have passed the threshold
betas_act = betaTrainAct(2:end); % remove the first weight which is the constant term

% sort interactors by decreasing weight and collect their names
betas_act(:,2) = 1:numel(betas_act); 
betas_act = sortrows(betas_act,'descend'); 
topPredictorsActivity = intList(InteractorIndAct(betas_act(1:nInteractorsToOutput,2)));


%% intensity
InteractorIndInt = find(sum(M)>=minInteractionsI); % retain only interactors that have passed the threshold
betas_int = betaTrainInt(2:end); % remove the first weight which is the constant term

% sort interactors by decreasing weight and collect their names
betas_int(:,2) = 1:numel(betas_int);
betas_int = sortrows(betas_int,'descend'); 
topPredictorsIntensity = intList(InteractorIndInt(betas_int(1:nInteractorsToOutput,2)));

end