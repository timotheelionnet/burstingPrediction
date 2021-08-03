function [M2,intList2,isInteractorInBothLists,tfList,lonelyTFs,discardedInts] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData,...
    minInteractionsPerTF,minInteractionsPerInteractor)
% remove interactors that do not make interactions with the characterized
% and non-characterized group

%% check that each interactor makes at least one interaction with a
% characterized TF and at least one interaction with a non-characterized
% TF.
isChar = ismember(tfList,burstingData.TFnames);

isNonChar = ~isChar;

isInteractorInBothLists = sum(M(isChar,:))>0 & sum(M(isNonChar,:))>0;

M2 = M(:,isInteractorInBothLists);
intList2 = intList(isInteractorInBothLists);

%% remove interactors that have fewer than minInteractionsPerInteractor
hasInteractorEnoughInteractions = sum(M2) >= minInteractionsPerInteractor;
discardedInts = intList2(~hasInteractorEnoughInteractions);
M2 = M2(:,hasInteractorEnoughInteractions);
intList2 = intList2(hasInteractorEnoughInteractions);

%% remove TFs that do not have any interactors (or fewer interactors than minInteractionsPerTF)
hasTFEnoughInteractions = sum(M2,2) >= minInteractionsPerTF;
lonelyTFs = tfList(~hasTFEnoughInteractions); 
M2 = M2(hasTFEnoughInteractions,:);
tfList = tfList(hasTFEnoughInteractions,:);




end