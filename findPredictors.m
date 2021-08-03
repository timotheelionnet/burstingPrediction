function [topPredictors_activity, topPredictors_intensity] ...
    = findPredictors(betaTrainAct, betaTrainInt, M, intList, A, I)
summation = sum(M);

betas_act = betaTrainAct;
betas_act(1,:) = [];
betas_act(:,2) = 1:numel(betas_act);
betas_act = sortrows(betas_act,'descend');
actSum = find(summation>=A);
topPredictors_activity = intList(actSum(betas_act(1:100,2)));

betas_int = betaTrainInt;
betas_int(1,:) = [];
betas_int(:,2) = 1:numel(betas_int);
betas_int = sortrows(betas_int,'descend');
intSum = find(summation>=I);
topPredictors_intensity = intList(intSum(betas_int(1:100,2)));

end