function [M2,intList2] = cutoffAndNormalizeInteractionsMatrix(M,intList,minInteractions)

% remove interactors that do not have at least minInteractions unique interactions
M2 = M(:,sum(M>0)>=minInteractions);
intList2 = intList(sum(M>0)>=minInteractions);

% z-score transform M2
M2 = zscore(zscore(M2')');

