%using the predicted AF and Int scores for all TFs, assign each TF to a
%cluster using the kmeans defined centroids and check for sig differences
%between clusters in protein features
figOpt = 1;
%% Initialize necessary data from BioGRID and LCR table
C = [1.62, 1.05; 1.02, 1.37; 2.7, 1.4; 0.99, 0.99; 1.6, 1.56];
C_names = {'AF','Int','High multi','Non','Low multi'};
if exist('RollingScores','var') == 0
    biogridScript; close all
end
if exist('allTF_lcr_table','var') == 0
    build_protein_features_table; close all
end
%% Find nearest centroid using euclidean distance
clusterID = zeros(length(tfList2),1);
for i = 1:length(tfList2)
    clusterID(i) = find(pdist2(RollingScores{5,2}(i,1:2),C) == min(pdist2(RollingScores{5,2}(i,1:2),C)));
end
%% Test for significant differences between clusters in protein features
h_AFvsX = zeros(2,width(allTF_lcr_table));
p_AFvsX = zeros(2,width(allTF_lcr_table));
for m = 1:length(C_names)
    for j = 2:width(allTF_lcr_table)
        [h_AFvsX(m,j),p_AFvsX(m,j)] = ttest2(table2array(allTF_lcr_table(find(clusterID == 1),j)),...
            table2array(allTF_lcr_table(find(clusterID == m),j)),'alpha',0.05/(width(allTF_lcr_table)-1));
    end
end
h_IntVsX = zeros(2,width(allTF_lcr_table));
p_IntVsX = zeros(2,width(allTF_lcr_table));
for m = 1:length(C_names)
    for j = 2:width(allTF_lcr_table)
        [h_IntVsX(m,j),p_IntVsX(m,j)] = ttest2(table2array(allTF_lcr_table(find(clusterID == 2),j)),...
            table2array(allTF_lcr_table(find(clusterID == m),j)),'alpha',0.05/(width(allTF_lcr_table)-1));
    end
end
h_HMvsX = zeros(2,width(allTF_lcr_table));
p_HMvsX = zeros(2,width(allTF_lcr_table));
for m = 1:numel(C_names)
    for j = 2:width(allTF_lcr_table)
        [h_HMvsX(m,j),p_HMvsX(m,j)] = ttest2(table2array(allTF_lcr_table(find(clusterID == 3),j)),...
            table2array(allTF_lcr_table(find(clusterID == m),j)),'alpha',0.05/(width(allTF_lcr_table)-1));
    end
end
h_NonVsX = zeros(2,width(allTF_lcr_table));
p_NonVsX = zeros(2,width(allTF_lcr_table));
for m = 1:numel(C_names)
    for j = 2:width(allTF_lcr_table)
        [h_NonVsX(m,j),p_NonVsX(m,j)] = ttest2(table2array(allTF_lcr_table(find(clusterID == 4),j)),...
            table2array(allTF_lcr_table(find(clusterID == m),j)),'alpha',0.05/(width(allTF_lcr_table)-1));
    end
end
h_LMvsX = zeros(2,width(allTF_lcr_table));
p_LMvsX = zeros(2,width(allTF_lcr_table));
for m = 1:numel(C_names)
    for j = 2:width(allTF_lcr_table)
        [h_LMvsX(m,j),p_LMvsX(m,j)] = ttest2(table2array(allTF_lcr_table(find(clusterID == 5),j)),...
            table2array(allTF_lcr_table(find(clusterID == m),j)),'alpha',0.05/(width(allTF_lcr_table)-1));
    end
end
%% Create bloxplots for features with sig differences, if plot flag is 1
if figOpt == 1
    for sigNum = 1:width(allTF_lcr_table)
        temp = sum(vertcat(h_AFvsX,h_IntVsX,h_HMvsX, h_NonVsX, h_LMvsX))>0;
        if temp(sigNum) == 1
            figure; boxplot(table2array(allTF_lcr_table(:,sigNum)),clusterID,'Notch','on')
            ylabel(allTF_lcr_table.Properties.VariableNames(sigNum) + " per TF")
            set(gca, 'xtick', 1:5); set(gca,'xticklabel',C_names)
        end
    end
end