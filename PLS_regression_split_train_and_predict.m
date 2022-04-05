figOpt = 0;
tic
itnum = 100;
A=240;
I=30;
InteractionIndexAct = find(sum(M3)>=A);
InteractionIndexInt = find(sum(M3)>=I);
M3act = M3norm(:,InteractionIndexAct);
M3int = M3norm(:,InteractionIndexInt);
%[m, n] = size(M3norm);
%for q = 1:100
%MscramAct = M3act(randperm(m),:);
%MscramInt = M3int(randperm(m),:);
trainingFraction = 10:10:80;
PLScorr_af = zeros(itnum,numel(trainingFraction));
PLScorr_int = zeros(itnum,numel(trainingFraction));
for i = 1:itnum
    for j = 1:numel(trainingFraction)
        nTraining = round(round(height(burstingData2)) * trainingFraction(j)/100);
        nValid = height(burstingData2) - nTraining;
        trainingIndex = randperm(round(height(burstingData2)),nTraining);
        validIndex = setdiff(1:height(burstingData2),trainingIndex);
        %% Train PLS model for AF effect
        [XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
            = plsregress(M3act(charTFind(trainingIndex),:),burstingData2.activity(trainingIndex),4);
        W0act = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
        %% Train PLS model for Int effect
        [XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
            = plsregress(M3int(charTFind(trainingIndex),:),burstingData2.intensity(trainingIndex),6);
        W0int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
        %% Predict AF scores for validation set
        predValidScoreAF = [ones(nValid,1), M3act(charTFind(validIndex),:)] * betaTrainAct;
        [ra,pa] = corr(burstingData2.activity(validIndex),predValidScoreAF,'type','Pearson');
        predValidTraincorrAct = ra;

        if figOpt == 1
            figure; scatter(burstingData2.activity(validIndex,1),predValidScoreAF); hold on;
            labelpoints(burstingData2.activity(validIndex,1),predValidScoreAF,burstingData2.TFnames(validIndex),'N',0.1)
            xlabel('True Activity Effect')
            ylabel('Predicted Activity Effect')
            title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
        end
        %% Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), M3int(charTFind(validIndex),:)] * betaTrainInt;
        [ri,pi] = corr(burstingData2.intensity(validIndex),predValidScoreInt,"type","Pearson");
        predValidTraincorrInt = ri;

        if figOpt == 1
            figure; scatter(burstingData2.intensity(validIndex,1),predValidScoreInt); hold on;
            labelpoints(burstingData2.intensity(validIndex,1),predValidScoreInt,burstingData2.TFnames(validIndex),'N',0.1)
            xlabel('True Intensity Effect')
            ylabel('Predicted Intensity Effect')
            title(['Predicted vs True Values - Intensity corr = ' num2str(predValidTraincorrInt)])
        end
        %%
        predTrainScoreact = [ones(nTraining,1), M3act(trainingIndex,:)] * betaTrainAct;
        predTrainTraincorrAct = corr(burstingData2.activity(trainingIndex),predTrainScoreact);
        predTrainScoreint = [ones(nTraining,1), M3int(trainingIndex,:)] * betaTrainInt;
        predTrainTraincorrInt = corr(burstingData2.intensity(trainingIndex),predTrainScoreint);
        %%
        PLScorr_af(i,j) = predValidTraincorrAct;
        PLScorr_int(i,j) = predValidTraincorrInt;

    end
    %end
    %PLScorr_af(q,:) = mean(PLScorr_af);
    %PLScorr_int(q,:) =  mean(PLScorr_int);
end
toc
figure;
errorbar(trainingFraction,mean(PLScorr_af),std(PLScorr_af)/sqrt(length(PLScorr_af)),'DisplayName','Active Fraction'); hold on
errorbar(trainingFraction,mean(PLScorr_int),std(PLScorr_int)/sqrt(length(PLScorr_int)),'DisplayName','Intensity');
xlabel('% of data used for training');
ylabel('mean PCC');
title("PLS");
legend
pbaspect([1 1 1])
box off