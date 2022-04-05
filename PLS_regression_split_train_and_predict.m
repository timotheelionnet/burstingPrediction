function [PLScorrAF,PLScorrInt,W0AF,W0Int] = PLS_regression_split_train_and_predict(...
    burstingData2,trainingFraction,M3norm,M3,charTFind,...
    minInteractionsA,minInteractionsI,nCompA,nCompI,itnum,figOpt)

tic

%% prune interactors
% select columns for interactors that have more than the threshold of
% required total interactions.
InteractionIndexAct = find(sum(M3)>=minInteractionsA);
InteractionIndexInt = find(sum(M3)>=minInteractionsI);

M3act = M3norm(:,InteractionIndexAct);
M3int = M3norm(:,InteractionIndexInt);

%%

%initializing the matrices holding correlation scores
PLScorrAF = zeros(itnum,numel(trainingFraction));
PLScorrInt = zeros(itnum,numel(trainingFraction));

for j = 1:numel(trainingFraction)
    disp(['training fraction ',num2str(trainingFraction(j)),'...']);
    for i = 1:itnum    
        % set the indices for training and validation sets
        nTraining = round(height(burstingData2) * trainingFraction(j)/100);
        nValid = height(burstingData2) - nTraining;
        trainingIndex = randperm(height(burstingData2),nTraining);
        validIndex = setdiff(1:height(burstingData2),trainingIndex);
        
        % Train PLS model for AF effect
        [~,~,~,~,betaTrainAct,~,~,statstrainact] = plsregress(...
            M3act(charTFind(trainingIndex),:),burstingData2.activity(trainingIndex),4);
        
        % compute
        W0AF = statstrainact.W ./ sqrt(sum(statstrainact.W.^2,1));
        
        % Train PLS model for Int effect
        [~,~,~,~,betaTrainInt,~,~,statstrainint] = plsregress(...
            M3int(charTFind(trainingIndex),:),burstingData2.intensity(trainingIndex),6);
        
        W0Int = statstrainint.W ./ sqrt(sum(statstrainint.W.^2,1));
        
        % Predict AF scores for validation set
        predValidScoreAF = [ones(nValid,1), M3act(charTFind(validIndex),:)] * betaTrainAct;
        
        % compute correlation for active fraction
        predValidTraincorrAct = corr(burstingData2.activity(validIndex),predValidScoreAF,...
            'type','Pearson');

%         if figOpt == 1
%             figure; hold on;
%             scatter(burstingData2.activity(validIndex,1),predValidScoreAF); 
%             labelpoints(burstingData2.activity(validIndex,1),...
%                 predValidScoreAF,burstingData2.TFnames(validIndex),'N',0.1);
%             xlabel('True Activity Effect')
%             ylabel('Predicted Activity Effect')
%             title(['Predicted vs True Values - Activity corr = ' num2str(predValidTraincorrAct)])
%         end
        
        % Predict Int scores for validation set
        predValidScoreInt = [ones(nValid,1), M3int(charTFind(validIndex),:)] * betaTrainInt;
        
        % compute correlation for intensity
        predValidTraincorrInt = corr(burstingData2.intensity(validIndex),predValidScoreInt,...
            "type","Pearson");
        
        % plot
%         if figOpt == 1
%             figure; hold on;  
%             scatter(burstingData2.intensity(validIndex,1),predValidScoreInt); 
%             labelpoints(burstingData2.intensity(validIndex,1),...
%                 predValidScoreInt,burstingData2.TFnames(validIndex),'N',0.1);
%             xlabel('True Intensity Effect')
%             ylabel('Predicted Intensity Effect')
%             title(['Predicted vs True Values - Intensity corr = ' num2str(predValidTraincorrInt)])
%         end
        
        %
        predTrainScoreact = [ones(nTraining,1), M3act(trainingIndex,:)] * betaTrainAct;
        predTrainTraincorrAct = corr(burstingData2.activity(trainingIndex),predTrainScoreact);
        predTrainScoreint = [ones(nTraining,1), M3int(trainingIndex,:)] * betaTrainInt;
        predTrainTraincorrInt = corr(burstingData2.intensity(trainingIndex),predTrainScoreint);
        
        % store correlation into matrix
        PLScorrAF(i,j) = predValidTraincorrAct;
        PLScorrInt(i,j) = predValidTraincorrInt;

    end
end
toc

%% plot correlations vs training fraction
figure('Name','Correlations between model prediction and measured validation set'); 
hold on
errorbar(trainingFraction,mean(PLScorrAF),...
    std(PLScorrAF)/sqrt(length(PLScorrAF)),...
    'DisplayName','Active Fraction'); 

errorbar(trainingFraction,mean(PLScorrInt),...
    std(PLScorrInt)/sqrt(length(PLScorrInt)),...
    'DisplayName','Intensity');

xlabel('% of data used for training');
ylabel('mean PCC');
title("PLS");
legend
pbaspect([1 1 1])
box off