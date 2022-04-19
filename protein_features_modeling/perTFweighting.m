if exist('RollingScores','var') == 0
    biogridScript; close all
end
if exist('lcr_table','var') == 0
    build_protein_features_table; close all
end
figOpt = 2;
%%
rp_act = zeros(width(allTF_lcr_table),2)';
rp_int = zeros(width(allTF_lcr_table),2)';
for i = 2:width(allTF_lcr_table)
    [rp_act(1,i),rp_act(2,i)] = corr(burstingData2.activity,table2array(allTF_lcr_table(charTFind,i)));
end
rp_act(3,:) = abs(rp_act(1,:))>0.1;

for i = 2:width(allTF_lcr_table)
    [rp_int(1,i),rp_int(2,i)] = corr(burstingData2.intensity,table2array(allTF_lcr_table(charTFind,i)));
end
rp_int(3,:) = abs(rp_int(1,:))>0.2;
%%
parameters = table2array(allTF_lcr_table(:,2:end));
parametersNorm = normalize(normalize(parameters,2));
parametersNorm(isnan(parametersNorm)) = 0;
SubAct = rp_act(3,2:end)>0;
parametersSubAct = parametersNorm(:,SubAct');
SubInt = rp_int(3,2:end)>0;
parametersSubInt = parametersNorm(:,SubInt');
%%
InteractionIndexAct = find(sum(M3)>=253);
InteractionIndexInt = find(sum(M3)>=40);
M3act = M3norm(charTFind,InteractionIndexAct);
M3int = M3norm(charTFind,InteractionIndexInt);
%%
itxPerTF = sum(M3(charTFind,:),2);
weight_act = 1.5;
weight_int = 4;
%%
ncompInteraction = [4,5];
ncompLCR = [1,1];
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
        trainingIndex = randperm(height(burstingData2),nTraining);
        validIndex = setdiff(1:height(burstingData2),trainingIndex);
        nValid = length(validIndex);
        predValidScoreAct = zeros(nValid,1);
        predValidScoreInt = zeros(nValid,1);
        %% Train PLS model for AF effect using interaction data
        [~,~,~,~,betaTrainActItx,~,~,~] ...
            = plsregress(M3act(trainingIndex,:),burstingData2.activity(trainingIndex),ncompInteraction(1));
        %% Train PLS model for Int effect using interaction data
        [~,~,~,~,betaTrainIntItx,~,~,~] ...
            = plsregress(M3int(trainingIndex,:),burstingData2.intensity(trainingIndex),ncompInteraction(2));
        %% Predict AF scores for validation set using interaction data
        predScores_act_interaction = [ones(nValid,1), M3act(validIndex,:)] * betaTrainActItx;
        %% Predict Int scores for validation set using interaction data
        predScores_int_interaction = [ones(nValid,1), M3int(validIndex,:)] * betaTrainIntItx;
        %% Train PLS model for AF effect using LCR data
        [~,~,~,~,betaTrainActLCR,~,~,~] ...
            = plsregress(parametersSubAct(trainingIndex,:),burstingData2.activity(trainingIndex),ncompLCR(1));
        %% Train PLS model for Int effect using LCR data
        [~,~,~,~,betaTrainIntLCR,~,~,~] ...
            = plsregress(parametersSubInt(trainingIndex,:),burstingData2.intensity(trainingIndex),ncompLCR(2));
        %% Predict AF scores for validation set using LCR data
        predScores_act_LCR = [ones(nValid,1), parametersSubAct(validIndex,:)] * betaTrainActLCR;
        %% Predict Int scores for validation set using LCR data
        predScores_int_LCR = [ones(nValid,1), parametersSubInt(validIndex,:)] * betaTrainIntLCR;
        %%
        for z = 1:nValid
            predValidScoreAct(z) = predScores_act_interaction(z)*((predScores_act_interaction(z) + weightInteractions(itxPerTF(z),...
                predScores_act_interaction(z),weight_act,median(itxPerTF)))/...
                (1+predScores_act_interaction(z) + weightInteractions(itxPerTF(z),predScores_act_interaction(z),weight_act,median(itxPerTF))))+...
                predScores_act_LCR(z)*(1/(1+predScores_act_LCR(z) + weightInteractions(itxPerTF(z),predScores_act_LCR(z),weight_act,median(itxPerTF))));
        end

        for z = 1:nValid
            predValidScoreInt(z) = predScores_int_interaction(z)*((predScores_int_interaction(z) + weightInteractions(itxPerTF(z),...
                predScores_int_interaction(z),weight_int,median(itxPerTF)))/...
                (1+predScores_int_interaction(z) + weightInteractions(itxPerTF(z),predScores_int_interaction(z),weight_int,median(itxPerTF))))+...
                predScores_act_LCR(z)*(1/(1+predScores_int_LCR(z) + weightInteractions(itxPerTF(z),predScores_int_LCR(z),weight_int,median(itxPerTF))));
        end

        predValidTraincorrAct = corr(burstingData2.activity(validIndex),predValidScoreAct);
        predValidTrainErrAct = immse(burstingData2.activity(validIndex), predValidScoreAct);
        predValidTraincorrInt = corr(burstingData2.intensity(validIndex),predValidScoreInt);
        predValidTrainErrInt = immse(burstingData2.intensity(validIndex), predValidScoreInt);
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