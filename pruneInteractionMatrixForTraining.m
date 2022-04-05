function [M2,intList2,isInteractorInBothLists,tfList2,lonelyTFs,discardedInts] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData,...
    minInteractionsPerTF,minInteractionsPerInteractor)
% From an interaction matrix M, this function generates a pruned interaction matrix M2, where 
% TFs and interactors that do not make enough total interactions have been removed. 
% Interactors that do not interact both with the experimentally characterized 
% and the non experimentally-characterized groups of of TFs are also removed from M. 

%%%%% INPUT
% M is the interaction matrix where the rows are TFs and the columns their interactors 
        % the size of M is <number of elements in tfList> x <number of elements in intList>
        % the genes in rows and columns of M follow the same order as in tfList and
        % intList
% intList: list of gene symbols for the columns of M (interactors)
% tfList: list of gene symbols for the rows of M (TFs)
% burstingData: MATLAB table with columns TFnames, activity, intensity (CRISPRburst results)  
        % the only column used here is 'TFnames'
% minInteractionsPerTF: % minimum number of interactions for a TF to be
        % kept in the pruned interaction matrix
% minInteractionsPerInteractor: % minimum number of interactions for an interactor to be
        % kept in the pruned interaction matrix

%%%%% OUTPUT
% M2: is the pruned interaction matrix where the rows are TFs and the columns their interactors
        % the size of M2 is <number of elements in tfList2> x <number of elements in intList2>
        % the genes in rows and columns of M2 follow the same order as in tfList2 and
        % intList2
% intList2: list of gene symbols for the columns of M2 (pruned interactors)
% isInteractorInBothLists: indices (relative to gene symbols in intList) of
        % the interactors (columns in M) that make interactions with: 
        % - at least one experimentally characterized TF, and
        % - at least one TF that is not experimentally characterized
% tfList: list of gene symbols for the rows of M2 (pruned TFs)
% lonelyTFs: list of of gene symbols for the TFs that were eventually 
        % discarded because they did not have more than
        % minInteractionsPerTF total interactions.
% discardedInts: list of gene symbols for the interactors that 
        % were interacting with at least one experimentally characterized TF
        % and at least one TF that is not experimentally characterized, 
        % BUT that were eventually discarded because they did not have more than
        % minInteractionsPerInteractor total interactions.

%% check that each interactor makes at least one interaction with a
% characterized TF and at least one interaction with a non-characterized
% TF.
isChar = ismember(tfList,burstingData.TFnames); % whether a TF (row in M) has been experimentally characterized

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
tfList2 = tfList(hasTFEnoughInteractions,:);

end