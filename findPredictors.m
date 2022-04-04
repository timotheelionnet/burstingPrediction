function [topPredictors_activity, topPredictors_intensity] ...
    = findPredictors(betaTrainAct, betaTrainInt, M, intList, A, I)

InteractionIndexStoreAct = find(sum(M)>=A);
betas_act = betaTrainAct(2:end);
betas_act(:,2) = 1:numel(betas_act);
betas_act = sortrows(betas_act,'descend');
topPredictors_activity = intList(InteractionIndexStoreAct(betas_act(1:100,2)));

InteractionIndexStoreInt = find(sum(M)>=I);
betas_int = betaTrainInt(2:end);
betas_int(:,2) = 1:numel(betas_int);
betas_int = sortrows(betas_int,'descend');
topPredictors_intensity = intList(InteractionIndexStoreInt(betas_int(1:100,2)));

end