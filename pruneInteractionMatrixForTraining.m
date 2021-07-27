function [M2,intList2,isInteractorInBothLists] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData)

isChar = ismember(tfList,burstingData.TFnames);

isNonChar = ~isChar;

isInteractorInBothLists = sum(M(isChar,:))>0 & sum(M(isNonChar,:))>0;

M2 = M(:,isInteractorInBothLists);
intList2 = intList(isInteractorInBothLists);



end