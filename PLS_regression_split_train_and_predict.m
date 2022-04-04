tic
itnum = 1000;
A=240;
I=30;
InteractionIndexAct = find(sum(M3)>=A);
InteractionIndexInt = find(sum(M3)>=I);
M3act = M3norm(:,InteractionIndexAct);
M3int = M3norm(:,InteractionIndexInt);
[m, n] = size(M3norm);
trainingFraction = 10:10:80;
for q = 1:100
MscramAct = M3act(randperm(m),:);
MscramInt = M3int(randperm(m),:);
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
        %figure; scatter(burstingData2.activity(validIndex,1),predValidScoreAF); hold on;
        %labelpoints(burstingData2.activity(validIndex,1),predValidScoreAF,charTFs(validIndex),'N',0.1)
        %xlabel('True Activity Effect')
        %ylabel('Predicted Activity Effect')
        [ra,pa,rla,ra] = corrcoef(burstingData2.activity(validIndex),predValidScoreAF);
        predValidTraincorrAct = ra;
        %title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
        %% Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), M3int(charTFind(validIndex),:)] * betaTrainInt;
        %figure; scatter(burstingData2.intensity(validIndex,1),predValidScoreInt); hold on;
        %labelpoints(burstingData2.intensity(validIndex,1),predValidScoreInt,charTFs(validIndex),'N',0.1)
        %xlabel('True Intensity Effect')
        %ylabel('Predicted Intensity Effect')
        [ri,pi,rli,ri] = corrcoef(burstingData2.intensity(validIndex),predValidScoreInt);
        predValidTraincorrInt = ri;
        %title(['Predicted vs True Values - Intensity corr = ' num2str(predValidTraincorrInt)])
        %%
        %predTrainScoreact = [ones(nTraining,1), M3act(charTFind,InteractionIndexAct)] * betaTrainAct;
        %predTrainTraincorrAct = corr(burstingData2.activity,predTrainScoreact);
        %predTrainScoreint = [ones(nTraining,1), M3int(charTFind,InteractionIndexInt)] * betaTrainInt;
        %predTrainTraincorrInt = corr(burstingData2.intensity,predTrainScoreint);
        %%
        PLScorr_af(i,j) = predValidTraincorrAct(1,2);
        PLScorr_int(i,j) = predValidTraincorrInt(1,2);

    end
end
PLScorr_af(q,:) = mean(PLScorr_af);
PLScorr_int(q,:) =  mean(PLScorr_int);
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