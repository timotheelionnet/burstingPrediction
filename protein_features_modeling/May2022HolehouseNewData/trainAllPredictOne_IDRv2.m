%% Run biogridScipt if needed
if exist('IDRfeatures','var') == 0
    PLS_regression_split_train_and_predict_IDRv2; close all
end
%%
figOpt = 1;
tic
PLSerror_af = [];
PLSerror_int = [];
predEff = [];
for j = 1:height(IDRfeatures)
    trainingIndex = setdiff(1:height(IDRfeatures),j);
    validIndex = j;
    nTraining = numel(trainingIndex);
    nValid = numel(validIndex);
    %% Train PLS model for AF effect
    [XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
        = plsregress(table2array(parametersSubAct(trainingIndex,:)),IDRfeatures.activity(trainingIndex),3);
        W0act = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
    %% Train PLS model for Int effect
    [XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
       = plsregress(table2array(parametersSubInt(trainingIndex,:)),IDRfeatures.intensity (trainingIndex),3);
        W0int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
    %% Predict AF scores for validation set
    predValidScoreAF = [ones(nValid,1), table2array(parametersSubAct(validIndex,:))] * betaTrainAct;
    predValidTrainErrorAct = abs(predValidScoreAF - IDRfeatures.activity(validIndex));
    predValidPercentErrorAct = abs(predValidScoreAF - IDRfeatures.activity(validIndex)) / IDRfeatures.activity(validIndex)*100;

    %% Predict Int scores for validation set
    predValidScoreInt = [ones(nValid,1), table2array(parametersSubInt(validIndex,:))] * betaTrainInt;
    predValidTrainErrorInt = abs(predValidScoreInt - IDRfeatures.intensity(validIndex));
    predValidPercentErrorInt = abs(predValidScoreInt - IDRfeatures.intensity(validIndex)) / IDRfeatures.intensity(validIndex)*100;

    %%
    predEff(j,1:2) = [predValidScoreAF,predValidScoreInt];
    PLSerror_af(j) = predValidTrainErrorAct;
    PLSerror_int(j) = predValidTrainErrorInt;
    PLSperEr_af(j) = predValidPercentErrorAct;
    PLSperEr_int(j) = predValidPercentErrorInt;
end
if figOpt == 1
toc
PLSerror_af = PLSerror_af';
PLSerror_af(:,2) = 1:length(PLSerror_af);
PLSerror_af = sortrows(PLSerror_af,1,'descend');
PLSerror_int = PLSerror_int';
PLSerror_int(:,2) = 1:length(PLSerror_int);
PLSerror_int = sortrows(PLSerror_int,1,'descend');
%%
figure; bar(PLSerror_af(:,1),'DisplayName','AF');
xlabel('TFs');
ylabel('abs error');
legend
set(gca,'XTickLabel',tfList2(charTFind(PLSerror_af(:,2))))
set(gca,'xtick',1:length(charTFind))
set(gca,'TickLength',[0 0])
xlim([0.5,73.5])
box off
figure; bar(PLSerror_int(:,1),'DisplayName','Int');
xlabel('TFs');
ylabel('abs error');
legend
set(gca,'XTickLabel',tfList2(charTFind(PLSerror_int(:,2))))
set(gca,'xtick',1:length(charTFind))
set(gca,'TickLength',[0 0])
xlim([0.5,73.5])
box off
end