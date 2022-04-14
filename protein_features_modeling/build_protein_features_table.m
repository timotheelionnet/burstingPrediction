%% Run biogridScipt if needed
if exist('RollingScores','var') == 0
    biogridScript; close all
end
%% Load LCR data
allTFproteinIDFileName = 'allTF UniProtKB IDs.txt';
allTFprotein_names  = readcell(allTFproteinIDFileName);
lcr_table_FileName = 'LCR_table.txt';
lcr_table  = readtable(lcr_table_FileName);
%% Create LCR table for all TF list
allTF_lcr_table(:,1) = table(tfList2);
allTF_lcr_table(:,2:23) = {0};
allTF_lcr_table.Properties.VariableNames = {'TFname','D','K','I','Y','G','R','M','E',...
    'L','W','P','F','H','T','N','A','C','V','S','Q','LCR cum length','LCR count'};
for i = 1:numel(allTFprotein_names)
    if strmatch(allTFprotein_names(i), lcr_table.Var1) ~= 0
        allTF_lcr_table(i,2:21) = array2table(mean(table2array(lcr_table(strmatch(string(allTFprotein_names(i))+'_', lcr_table.Var1),2:21)),1));
        allTF_lcr_table(i,22) = array2table(sum(table2array(lcr_table(strmatch(string(allTFprotein_names(i))+'_', lcr_table.Var1),24))));
        allTF_lcr_table(i,23) = array2table(max(table2array(lcr_table(strmatch(string(allTFprotein_names(i))+'_', lcr_table.Var1),25))));
    else
        allTF_lcr_table(i,2:23) = {0};
    end
end
%% Find signficant correlations between predicted effects and protein features with FDR correction (alpha / # tests)
rp_act = zeros(width(allTF_lcr_table),2)';
rp_int = zeros(width(allTF_lcr_table),2)';
for i = 2:width(allTF_lcr_table)
    [rp_act(1,i),rp_act(2,i)] = corr(RollingScores{5,2}(:,1),table2array(allTF_lcr_table(:,i)));
end
rp_act(3,:) = rp_act(2,:)<0.05/(length(rp_act')-1);
sigCorr_act = allTF_lcr_table.Properties.VariableNames(rp_act(3,:)>0);
for i = 2:width(allTF_lcr_table)
    [rp_int(1,i),rp_int(2,i)] = corr(RollingScores{5,2}(:,2),table2array(allTF_lcr_table(:,i)));
end
rp_int(3,:) = rp_int(2,:)<0.05/(length(rp_int')-1);
sigCorr_int = allTF_lcr_table.Properties.VariableNames(rp_int(3,:)>0);
%% Plot some significant correlations with amino acids
figure; scatter(allTF_lcr_table.Q,RollingScores{5,2}(:,1)); lsline;
ylim([0 3.5]); ylabel('Predicted Active Fraction Effect'); xlabel('Frequency of Glutamine per TF');

figure; scatter(allTF_lcr_table.G,RollingScores{5,2}(:,2)); lsline;
ylim([0.5 2]); ylabel('Predicted Intensity Effect'); xlabel('Frequency of Glycine per TF');

% Plot some significant correlations with LCR features
figure; scatter(allTF_lcr_table.("LCR count"),RollingScores{5,2}(:,2)); lsline;
ylim([0.5 2]); ylabel('Predicted Active Fraction Effect'); xlabel(allTF_lcr_table.Properties.VariableNames(23) + " per TF");

figure; scatter(allTF_lcr_table.("LCR cum length"),RollingScores{5,2}(:,2)); lsline;
ylim([0.5 2]); ylabel('Predicted Intensity Effect'); xlabel(allTF_lcr_table.Properties.VariableNames(22) + " per TF");