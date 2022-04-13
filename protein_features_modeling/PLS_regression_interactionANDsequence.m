if exist('RollingScores','var') == 0
    biogridScript; close all
end
if exist('lcr_table','var') == 0
    build_protein_features_table; close all
end
%%
parameters = table2array(allTF_lcr_table(:,2:end));
parametersNorm = normalize(normalize(parameters,2));
parametersNorm(isnan(parametersNorm)) = 0;
SubAct = rp_act(2,2:end)<0.05;
parametersSubAct = parametersNorm(:,SubAct');
SubInt = rp_int(2,2:end)<0.05;
parametersSubInt = parametersNorm(:,SubInt');
%%
figOpt = 2;
A=270;
I=40;
M3act = M3norm(:,sum(M3)>=A);
M3int = M3norm(:,sum(M3)>=I);
ActData = [M3act(charTFind,:),parametersSubAct(charTFind,:)];
IntData = [M3int(charTFind,:),parametersSubInt(charTFind,:)];

ncomp = [2,3];
itnum = 100;
trainingFraction = 10:10:80;

PLScorr_af = zeros(itnum,numel(trainingFraction));
PLSerr_af = zeros(itnum,numel(trainingFraction));
PLScorr_int = zeros(itnum,numel(trainingFraction));
PLSerr_int = zeros(itnum,numel(trainingFraction));
tic
%%
for i = 1:itnum
    for j = 1:numel(trainingFraction)
        nTraining = round(round(height(burstingData2)) * trainingFraction(j)/100); 
        %keepIndex = [4,6:73];
        trainingIndex = randperm(height(burstingData2),nTraining);
        %trainingIndex = trainingIndex(ismember(trainingIndex,keepIndex));
        validIndex = setdiff(1:height(burstingData2),trainingIndex);
        %validIndex = validIndex(ismember(validIndex,keepIndex));
        nValid = length(validIndex);
        %% Train PLS model for AF effect
        [XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
            = plsregress(ActData(trainingIndex,:),burstingData2.activity(trainingIndex),ncomp(1));
        W0act = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
        %% Train PLS model for Int effect
        [XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
            = plsregress(IntData(trainingIndex,:),burstingData2.intensity(trainingIndex),ncomp(2));
        W0int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
        %% Predict AF scores for validation set
        predValidScoreAct = [ones(nValid,1), ActData(validIndex,:)] * betaTrainAct;
        predValidTraincorrAct = corr(burstingData2.activity(validIndex),predValidScoreAct);
        predValidTrainErrAct = immse(burstingData2.activity(validIndex), predValidScoreAct);

        if figOpt == 1
            figure; scatter(burstingData2.activity(validIndex,1),predValidScoreAct); hold on;
            labelpoints(burstingData2.activity(validIndex,1),predValidScoreAct,burstingData2.TFnames(validIndex),'N',0.1)
            xlabel('True Activity Effect')
            ylabel('Predicted Activity Effect')
            title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
        end
        %% Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), IntData(validIndex,:)] * betaTrainInt;
        predValidTraincorrInt = corr(burstingData2.intensity(validIndex),predValidScoreInt);
        predValidTrainErrInt = immse(burstingData2.intensity(validIndex), predValidScoreInt);


        if figOpt == 1
            figure; scatter(burstingData2.intensity(validIndex,1),predValidScoreInt); hold on;
            labelpoints(burstingData2.intensity(validIndex,1),predValidScoreInt,burstingData2.TFnames(validIndex),'N',0.1)
            xlabel('True Intensity Effect')
            ylabel('Predicted Intensity Effect')
            title(['Predicted vs True Values - Intensity corr = ' num2str(predValidTraincorrInt)])
        end
        %% Collect metrics
        PLScorr_af(i,j) = predValidTraincorrAct;
        PLSerr_af(i,j) = mean(predValidTrainErrAct);
        PLScorr_int(i,j) = predValidTraincorrInt;
        PLSerr_int(i,j) = mean(predValidTrainErrInt);
    end
end
toc
if figOpt == 1 || 2
    figure;
    errorbar(trainingFraction,mean(PLScorr_af),std(PLScorr_af)/sqrt(length(PLScorr_af)),'DisplayName','Active Fraction'); hold on
    errorbar(trainingFraction,mean(PLScorr_int),std(PLScorr_int)/sqrt(length(PLScorr_int)),'DisplayName','Intensity');
    xlabel('trainingFraction');
    ylabel(['mean PCC over ' num2str(itnum) ' iterations']);
    title("PLS interaction+protein");
    legend(location = 'best'); pbaspect([1 1 1]); box off

    figure;
    errorbar(trainingFraction,mean(PLSerr_af),std(PLSerr_af)/sqrt(length(PLSerr_af)),'DisplayName','Active Fraction'); hold on
    errorbar(trainingFraction,mean(PLSerr_int),std(PLSerr_int)/sqrt(length(PLSerr_int)),'DisplayName','Intensity');
    xlabel('trainingFraction');
    ylabel(['MSE over ' num2str(itnum) ' iterations']);
    title("PLS interaction+protein");
    legend(location = 'best'); pbaspect([1 1 1]); box off
end