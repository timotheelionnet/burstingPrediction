figOpt = 2;
IDRdataFileName = 'norm_parameters - FullProtein.csv';
%IDRdataFileName = 'norm_parameters - allIDR.csv';
%IDRdataFileName = 'norm_parameters - 60IDR.csv';
IDRfeatures = readtable(IDRdataFileName);
%%
rp_act = zeros(width(IDRfeatures)-3,2)';
rp_int = zeros(width(IDRfeatures)-3,2)';

for i = 1:width(IDRfeatures)-3
    [rp_act(1,i),rp_act(2,i)] = corr(IDRfeatures.activity,table2array(IDRfeatures(:,i+3)));
end
for i = 1:width(IDRfeatures)-3
    [rp_int(1,i),rp_int(2,i)] = corr(IDRfeatures.intensity,table2array(IDRfeatures(:,i+3)));
end
rp_act(3,:) = rp_act(1,:)>0.05;
rp_int(3,:) = rp_int(1,:)>0.05;
%%
SubAct = logical([0,0,0,rp_act(3,:)>0]);
parametersSubAct = IDRfeatures(:,SubAct');
SubInt = logical([0,0,0,rp_int(3,:)>0]);
parametersSubInt = IDRfeatures(:,SubInt');

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
        nTraining = round(round(height(parametersSubAct)) * trainingFraction(j)/100);
        %keepIndex = [4,6:73];
        trainingIndex = randperm(height(parametersSubAct),nTraining);
        %trainingIndex = trainingIndex(ismember(trainingIndex,keepIndex));
        validIndex = setdiff(1:height(parametersSubAct),trainingIndex);
        %validIndex = validIndex(ismember(validIndex,keepIndex));
        nValid = length(validIndex);
        %% Train PLS model for AF effect
        [XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
            = plsregress(table2array(parametersSubAct(trainingIndex,:)),IDRfeatures.activity(trainingIndex),3);
        W0act = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
        %% Train PLS model for Int effect
        [XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
            = plsregress(table2array(parametersSubInt(trainingIndex,:)),IDRfeatures.intensity(trainingIndex),3);
        W0int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
        %% Predict AF scores for validation set
        predValidScoreAct = [ones(nValid,1), table2array(parametersSubAct(validIndex,:))] * betaTrainAct;
        predValidTraincorrAct = corr(IDRfeatures.activity(validIndex),predValidScoreAct);
        predValidTrainErrAct = immse(IDRfeatures.activity(validIndex), predValidScoreAct);
        predValidTrainMedErrAct = median(IDRfeatures.activity(validIndex)) - median(predValidScoreAct);

        if figOpt == 1
            figure; scatter(IDRfeatures.activity(validIndex),predValidScoreAct); hold on;
            labelpoints(IDRfeatures.activity(validIndex),predValidScoreAct,tfList2(charTFind(validIndex)),'N',0.1)
            xlabel('True Activity Effect')
            ylabel('Predicted Activity Effect')
            title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
        end
        %% Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), table2array(parametersSubInt(validIndex,:))] * betaTrainInt;
        predValidTraincorrInt = corr(IDRfeatures.intensity(validIndex),predValidScoreInt);
        predValidTrainErrInt = immse(IDRfeatures.intensity(validIndex), predValidScoreInt);
        predValidTrainMedErrInt = median(IDRfeatures.intensity(validIndex)) - median(predValidScoreInt);


        if figOpt == 1
            figure; scatter(IDRfeatures.intensity(validIndex),predValidScoreInt); hold on;
            labelpoints(IDRfeatures.intensity(validIndex),predValidScoreInt,tfList2(charTFind(validIndex)),'N',0.1)
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
    title("PLS protein features only");
    legend(location = 'best'); pbaspect([1 1 1]); box off

    figure;
    errorbar(trainingFraction,mean(PLSerr_afpro),std(PLSerr_afpro)/sqrt(length(PLSerr_afpro)),'DisplayName','Active Fraction'); hold on
    errorbar(trainingFraction,mean(PLSerr_intpro),std(PLSerr_intpro)/sqrt(length(PLSerr_intpro)),'DisplayName','Intensity');
    xlabel('% of data used for training');
    ylabel(['MSE over ' num2str(itnum) ' iterations']);
    title("PLS protein features only");
    legend(location = 'best'); pbaspect([1 1 1]); box off

    %figure;
    %errorbar(trainingFraction,mean(PLSmedianErr_afpro),std(PLSmedianErr_afpro)/sqrt(length(PLSmedianErr_afpro)),'DisplayName','Active Fraction'); hold on
    %errorbar(trainingFraction,mean(PLSmedianErr_intpro),std(PLSmedianErr_intpro)/sqrt(length(PLSmedianErr_intpro)),'DisplayName','Intensity');
    %xlabel('% of data used for training');
    %ylabel(['Error between median TF predictions','(true - pred)']);
    %title("repressors in validation only");
    %legend(location = 'best'); pbaspect([1 1 1]); box off
end