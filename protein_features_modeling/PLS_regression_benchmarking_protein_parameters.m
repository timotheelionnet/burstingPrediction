figOpt = 2;
parametersNorm = normalize(normalize(parameters,2));
parametersNorm(isnan(parametersNorm)) = 0;

rp_act = zeros(width(allTF_lcr_table),2)';
rp_int = zeros(width(allTF_lcr_table),2)';
for i = 2:width(allTF_lcr_table)
    [rp_act(1,i),rp_act(2,i)] = corr(RollingScores{5,2}(:,1),table2array(allTF_lcr_table(:,i)));
end
SubAct = rp_act(2,2:end)<0.05;
parametersSubAct = parametersNorm(charTFind,SubAct');

for i = 2:width(allTF_lcr_table)
    [rp_int(1,i),rp_int(2,i)] = corr(RollingScores{5,2}(:,2),table2array(allTF_lcr_table(:,i)));
end
SubInt = rp_int(2,2:end)<0.05;
parametersSubInt = parametersNorm(charTFind,SubInt');

ncomp = [3,3];
tic
itnum = 100;
trainingFraction = 10:10:80;
PLScorr_afpro = zeros(itnum,numel(trainingFraction));
PLSerr_afpro = zeros(itnum,numel(trainingFraction));
PLScorr_intpro = zeros(itnum,numel(trainingFraction));
PLSerr_intpro = zeros(itnum,numel(trainingFraction));
PLSmedianErr_afpro = zeros(itnum,numel(trainingFraction));
PLSmedianErr_intpro = zeros(itnum,numel(trainingFraction));
for i = 1:itnum
    for j = 1:numel(trainingFraction)
        nTraining = round(round(length(parametersSubAct)) * trainingFraction(j)/100);
        %keepIndex = [4,6:73];
        trainingIndex = randperm(length(parametersSubAct),nTraining);
        %trainingIndex = trainingIndex(ismember(trainingIndex,keepIndex));
        validIndex = setdiff(1:length(parametersSubAct),trainingIndex);
        %validIndex = validIndex(ismember(validIndex,keepIndex));
        nValid = length(validIndex);
        %% Train PLS model for AF effect
        [XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
            = plsregress(parametersSubAct(trainingIndex,:),burstingData2.activity(trainingIndex),ncomp(1));
        W0act = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
        %% Train PLS model for Int effect
        [XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
            = plsregress(parametersSubInt(trainingIndex,:),burstingData2.intensity (trainingIndex),ncomp(2));
        W0int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
        %% Predict AF scores for validation set
        predValidScoreAct = [ones(nValid,1), parametersSubAct(validIndex,:)] * betaTrainAct;
        predValidTraincorrAct = corr(burstingData2.activity(validIndex),predValidScoreAct);
        predValidTrainErrAct = immse(burstingData2.activity(validIndex), predValidScoreAct);
        predValidTrainMedErrAct = median(burstingData2.activity(validIndex)) - median(predValidScoreAct);

        if figOpt == 1
            figure; scatter(burstingData2.activity(validIndex,1),predValidScoreAct); hold on;
            labelpoints(burstingData2.activity(validIndex,1),predValidScoreAct,tfList2(charTFind(validIndex)),'N',0.1)
            xlabel('True Activity Effect')
            ylabel('Predicted Activity Effect')
            title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
        end
        %% Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), parametersSubInt(validIndex,:)] * betaTrainInt;
        predValidTraincorrInt = corr(burstingData2.intensity(validIndex),predValidScoreInt);
        predValidTrainErrInt = immse(burstingData2.intensity(validIndex), predValidScoreInt);
        predValidTrainMedErrInt = median(burstingData2.intensity(validIndex)) - median(predValidScoreInt);


        if figOpt == 1
            figure; scatter(burstingData2.intensity(validIndex,1),predValidScoreInt); hold on;
            labelpoints(burstingData2.intensity(validIndex,1),predValidScoreInt,tfList2(charTFind(validIndex)),'N',0.1)
            xlabel('True Intensity Effect')
            ylabel('Predicted Intensity Effect')
            title(['Predicted vs True Values - Intensity corr = ' num2str(predValidTraincorrInt)])
        end
        %% Collect metrics
        PLScorr_afpro(i,j) = predValidTraincorrAct;
        PLSerr_afpro(i,j) = mean(predValidTrainErrAct);
        PLScorr_intpro(i,j) = predValidTraincorrInt;
        PLSerr_intpro(i,j) = mean(predValidTrainErrInt);
        PLSmedianErr_afpro(i,j) = predValidTrainMedErrAct;
        PLSmedianErr_intpro(i,j) = predValidTrainMedErrInt;
    end
end
toc
if figOpt == 1 || 2
    figure;
    errorbar(trainingFraction,mean(PLScorr_afpro),std(PLScorr_afpro)/sqrt(length(PLScorr_afpro)),'DisplayName','Active Fraction'); hold on
    errorbar(trainingFraction,mean(PLScorr_intpro),std(PLScorr_intpro)/sqrt(length(PLScorr_intpro)),'DisplayName','Intensity');
    xlabel('% of data used for training');
    ylabel(['mean PCC over ' num2str(itnum) ' iterations']);
    title("repressors in validation only");
    legend(location = 'best'); pbaspect([1 1 1]); box off

    figure;
    errorbar(trainingFraction,mean(PLSerr_afpro),std(PLSerr_afpro)/sqrt(length(PLSerr_afpro)),'DisplayName','Active Fraction'); hold on
    errorbar(trainingFraction,mean(PLSerr_intpro),std(PLSerr_intpro)/sqrt(length(PLSerr_intpro)),'DisplayName','Intensity');
    xlabel('% of data used for training');
    ylabel(['MSE over ' num2str(itnum) ' iterations']);
    title("repressors in validation only");
    legend(location = 'best'); pbaspect([1 1 1]); box off

    %figure;
    %errorbar(trainingFraction,mean(PLSmedianErr_afpro),std(PLSmedianErr_afpro)/sqrt(length(PLSmedianErr_afpro)),'DisplayName','Active Fraction'); hold on
    %errorbar(trainingFraction,mean(PLSmedianErr_intpro),std(PLSmedianErr_intpro)/sqrt(length(PLSmedianErr_intpro)),'DisplayName','Intensity');
    %xlabel('% of data used for training');
    %ylabel(['Error between median TF predictions','(true - pred)']);
    %title("repressors in validation only");
    %legend(location = 'best'); pbaspect([1 1 1]); box off
end