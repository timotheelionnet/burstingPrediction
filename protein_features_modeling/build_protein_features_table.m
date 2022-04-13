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
    'L','W','P','F','H','T','N','A','C','V','S','Q','LCR_cum_length','LCR_count'};
for i = 1:numel(allTFprotein_names)
    if strmatch(allTFprotein_names(i), lcr_table.Var1) ~= 0
        allTF_lcr_table(i,2:21) = array2table(mean(table2array(lcr_table(strmatch(allTFprotein_names(i), lcr_table.Var1),2:21)),1));
        allTF_lcr_table(i,22) = array2table(sum(table2array(lcr_table(strmatch(allTFprotein_names(i), lcr_table.Var1),24))));
        allTF_lcr_table(i,23) = array2table(max(table2array(lcr_table(strmatch(allTFprotein_names(i), lcr_table.Var1),25))));
    end
end
%% Find signficant correlations between predicted effects and protein features with FDR correction
rp_act = zeros(width(allTF_lcr_table),2)';
rp_int = zeros(width(allTF_lcr_table),2)';
for i = 2:width(allTF_lcr_table)
    [rp_act(1,i),rp_act(2,i)] = corr(RollingScores{5,2}(:,1),table2array(allTF_lcr_table(:,i)));
end
rp_act(3,:) = mafdr(rp_act(2,:));
rp_act(4,:) = rp_act(2,:)<rp_act(3,:);
for i = 2:width(allTF_lcr_table)
    [rp_int(1,i),rp_int(2,i)] = corr(RollingScores{5,2}(:,2),table2array(allTF_lcr_table(:,i)));
end
rp_int(3,:) = mafdr(rp_int(2,:));
rp_int(4,:) = rp_int(2,:)<rp_int(3,:);
%% Plot signficant correlations after FDR correction
figure; scatter(table2array(allTF_lcr_table(:,3)),RollingScores{5,2}(:,1))
ylim([0 3.5]); ylabel('Predicted Active Fraction Effect'); xlabel('Frequency of Lysine per TF');
figure; scatter(table2array(allTF_lcr_table(:,22)),RollingScores{5,2}(:,1))
ylim([0 3.5]); ylabel('Predicted Active Fraction Effect'); xlabel('Cumulative length of LCR segments per TF');
figure; scatter(table2array(allTF_lcr_table(:,23)),RollingScores{5,2}(:,1))
ylim([0 3.5]); ylabel('Predicted Active Fraction Effect'); xlabel('Number of LCR segments per TF');
