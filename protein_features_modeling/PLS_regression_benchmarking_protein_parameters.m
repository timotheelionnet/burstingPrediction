figOpt = 2;
parametersSub = parameters(:,1:48);
parametersSub = normalize(parametersSub,"range");
ncomp = 2;
tic
itnum = 100;
trainingFraction = 10:10:80;
PLScorr_afpro = zeros(itnum,numel(trainingFraction));
PLSerr_afpro = zeros(itnum,numel(trainingFraction));
PLScorr_intpro = zeros(itnum,numel(trainingFraction));
PLSerr_intpro = zeros(itnum,numel(trainingFraction));
for i = 1:itnum
    for j = 1:numel(trainingFraction)
        nTraining = round(round(height(burstingDataParameters)) * trainingFraction(j)/100);
        nValid = height(burstingDataParameters) - nTraining;
        trainingIndex = randperm(round(height(burstingDataParameters)),nTraining);
        validIndex = setdiff(1:height(burstingDataParameters),trainingIndex);
        %% Train PLS model for AF effect
        [XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
            = plsregress(parametersSub(trainingIndex,:),burstingDataParameters.activity(trainingIndex),ncomp);
        W0act = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
        %% Train PLS model for Int effect
        [XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
            = plsregress(parametersSub(trainingIndex,:),burstingDataParameters.intensity (trainingIndex),ncomp);
        W0int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
        %% Predict AF scores for validation set
        predValidScoreAF = [ones(nValid,1), parametersSub(validIndex,:)] * betaTrainAct;
        predValidTraincorrAct = corr(burstingDataParameters.activity(validIndex),predValidScoreAF);
        predValidTrainErrAct = immse(burstingDataParameters.activity(validIndex), predValidScoreAF);

        if figOpt == 1
            figure; scatter(burstingDataParameters.activity(validIndex,1),predValidScoreAF); hold on;
            labelpoints(burstingDataParameters.activity(validIndex,1),predValidScoreAF,parametersNames(validIndex),'N',0.1)
            xlabel('True Activity Effect')
            ylabel('Predicted Activity Effect')
            title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
        end
        %% Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), parametersSub(validIndex,:)] * betaTrainInt;
        predValidTraincorrInt = corr(burstingDataParameters.intensity(validIndex),predValidScoreInt);
        predValidTrainErrInt = immse(burstingDataParameters.intensity(validIndex), predValidScoreInt);

        if figOpt == 1
            figure; scatter(burstingDataParameters.intensity(validIndex,1),predValidScoreInt); hold on;
            labelpoints(burstingDataParameters.intensity(validIndex,1),predValidScoreInt,parametersNames(validIndex),'N',0.1)
            xlabel('True Intensity Effect')
            ylabel('Predicted Intensity Effect')
            title(['Predicted vs True Values - Intensity corr = ' num2str(predValidTraincorrInt)])
        end
        PLScorr_afpro(i,j) = predValidTraincorrAct;
        PLSerr_afpro(i,j) = mean(predValidTrainErrAct);
        PLScorr_intpro(i,j) = predValidTraincorrInt;
        PLSerr_intpro(i,j) = mean(predValidTrainErrInt);
    end
end
toc
if figOpt == 1 || 2
    figure;
    errorbar(trainingFraction,mean(PLScorr_afpro),std(PLScorr_afpro)/sqrt(length(PLScorr_afpro)),'DisplayName','Active Fraction'); hold on
    errorbar(trainingFraction,mean(PLScorr_intpro),std(PLScorr_intpro)/sqrt(length(PLScorr_intpro)),'DisplayName','Intensity');
    xlabel('% of data used for training');
    ylabel(['mean PCC over ' num2str(itnum) ' iterations']);
    title("PLS protein");
    legend(location = 'best'); pbaspect([1 1 1]); box off
    figure;
    errorbar(trainingFraction,mean(PLSerr_afpro),std(PLSerr_afpro)/sqrt(length(PLSerr_afpro)),'DisplayName','Active Fraction'); hold on
    errorbar(trainingFraction,mean(PLSerr_intpro),std(PLSerr_intpro)/sqrt(length(PLSerr_intpro)),'DisplayName','Intensity');
    xlabel('% of data used for training');
    ylabel(['MSE over ' num2str(itnum) ' iterations']);
    title("PLS protein");
    legend(location = 'best'); pbaspect([1 1 1]); box off
end