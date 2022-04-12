figOpt = 0;

ransacOpt = 1;
if ransacOpt == 1
    effects = table2array(burstingData2(:,2:3));
    sampleSize = 2; % number of points to sample per trial
    maxDistance = 0.5; % max allowable distance for inliers
    fitLineFcn = @(effects) polyfit(effects(:,1),effects(:,2),1);
    evalLineFcn = @(model, effects) sum((effects(:,2)-polyval(model,effects(:,1))).^2,2); 
    ncomp = [2,2];
else
    ncomp = [4,5];
end

itnum = 100;
A=253;
I=40;
InteractionIndexAct = find(sum(M3)>=A);
InteractionIndexInt = find(sum(M3)>=I);
M3act = M3norm(:,InteractionIndexAct);
M3int = M3norm(:,InteractionIndexInt);
trainingFraction = 10:10:80;
PLScorr_af = zeros(itnum,numel(trainingFraction));
PLScorr_int = zeros(itnum,numel(trainingFraction));
PLSerr_af = zeros(itnum,numel(trainingFraction));
PLSerr_int = zeros(itnum,numel(trainingFraction));

tic
for i = 1:itnum
    for j = 1:numel(trainingFraction)
        if ransacOpt == 1
            [modelRANSAC, inlierIdx] = ransac(effects,fitLineFcn,evalLineFcn, ...
                sampleSize,maxDistance);
            outlierIdx = find(inlierIdx == 0);
        else
            inlierIdx = 1:73;
        end
        burstingData2sub = burstingData2(inlierIdx,:);
        nTraining = round(round(height(burstingData2sub)) * trainingFraction(j)/100);
        nValid = height(burstingData2sub) - nTraining;
        trainingIndex = randperm(round(height(burstingData2sub)),nTraining);
        validIndex = setdiff(1:height(burstingData2sub),trainingIndex);
        %% Train PLS model for AF effect
        [XLtrainact,yltrainact,XStrainact,YStrainact,betaTrainAct,PCTVARtrainact,splitMSEtrainact,statstrainact] ...
            = plsregress(M3act(charTFind(trainingIndex),:),burstingData2sub.activity(trainingIndex),ncomp(1));
        W0act = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
        %% Train PLS model for Int effect
        [XLtrainint,yltrainint,XStrainint,YStrainint,betaTrainInt,PCTVARtrainint,splitMSEtrainint,statstrainint] ...
            = plsregress(M3int(charTFind(trainingIndex),:),burstingData2sub.intensity(trainingIndex),ncomp(2));
        W0int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
        %% Predict AF scores for validation set
        predValidScoreAF = [ones(nValid,1), M3act(charTFind(validIndex),:)] * betaTrainAct;
        predValidTraincorrAct = corr(burstingData2sub.activity(validIndex),predValidScoreAF,'type','Pearson');
        predValidTrainErrAct = immse(burstingData2sub.activity(validIndex),predValidScoreAF);

        if figOpt == 1
            figure; scatter(burstingData2sub.activity(validIndex,1),predValidScoreAF); hold on;
            labelpoints(burstingData2sub.activity(validIndex,1),predValidScoreAF,burstingData2sub.TFnames(validIndex),'N',0.1)
            xlabel('True Activity Effect')
            ylabel('Predicted Activity Effect')
            title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
        end
        %% Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), M3int(charTFind(validIndex),:)] * betaTrainInt;
        predValidTraincorrInt = corr(burstingData2sub.intensity(validIndex),predValidScoreInt,"type","Pearson");
        predValidTrainErrInt = immse(burstingData2sub.intensity(validIndex),predValidScoreInt);

        if figOpt == 1
            figure; scatter(burstingData2sub.intensity(validIndex,1),predValidScoreInt); hold on;
            labelpoints(burstingData2sub.intensity(validIndex,1),predValidScoreInt,burstingData2sub.TFnames(validIndex),'N',0.1)
            xlabel('True Intensity Effect')
            ylabel('Predicted Intensity Effect')
            title(['Predicted vs True Values - Intensity corr = ' num2str(predValidTraincorrInt)])
        end
        %%
        predTrainScoreact = [ones(nTraining,1), M3act(trainingIndex,:)] * betaTrainAct;
        predTrainTraincorrAct = corr(burstingData2sub.activity(trainingIndex),predTrainScoreact);

        predTrainScoreint = [ones(nTraining,1), M3int(trainingIndex,:)] * betaTrainInt;
        predTrainTraincorrInt = corr(burstingData2sub.intensity(trainingIndex),predTrainScoreint);

        %%
        PLScorr_af(i,j) = predValidTraincorrAct;
        PLScorr_int(i,j) = predValidTraincorrInt;
        PLSerr_af(i,j) = predValidTrainErrAct;
        PLSerr_int(i,j) = predValidTrainErrInt;
        if ransacOpt == 1
        inlierCollect(1:73,i) = inlierIdx;
        end
    end
end
if ransacOpt == 1
outlierProb = 1-mean(inlierCollect,2);
outlierProb(:,2) = 1:length(outlierProb);
outlierProb = sortrows(outlierProb, "descend");
figure;
bar(outlierProb(:,1))
set(gca,'xticklabel',burstingData2.TFnames(outlierProb(:,2)))
ylabel("Outlier probability",'FontSize', 14);
set(gca,'xtick',1:height(burstingData2))
xtickangle(40); set(gca, 'ticklength', [0,0]);
box off
end
toc
figure;
errorbar(trainingFraction,mean(PLScorr_af),std(PLScorr_af)/sqrt(length(PLScorr_af)),'DisplayName','Active Fraction'); hold on
errorbar(trainingFraction,mean(PLScorr_int),std(PLScorr_int)/sqrt(length(PLScorr_int)),'DisplayName','Intensity');
xlabel('% of data used for training');
ylabel('PCC');
title("PLS");
legend
pbaspect([1 1 1])
box off

figure;
errorbar(trainingFraction,mean(PLSerr_af),std(PLSerr_af)/sqrt(length(PLSerr_af)),'DisplayName','Active Fraction'); hold on
errorbar(trainingFraction,mean(PLSerr_int),std(PLSerr_int)/sqrt(length(PLSerr_int)),'DisplayName','Intensity');
xlabel('% of data used for training');
ylabel('MSE');
title("PLS");
legend
pbaspect([1 1 1])
box off