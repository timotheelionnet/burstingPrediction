figOpt = 1;
IDRdataFileName = 'norm_parameters - FullProtein.csv';
%IDRdataFileName = 'norm_parameters - allIDR.csv';
%IDRdataFileName = 'norm_parameters - 60IDR.csv';
IDRfeatures = readtable(IDRdataFileName);
%%
C = [1.62, 1.05; 1.02, 1.37; 2.7, 1.4; 0.99, 0.99; 1.6, 1.56];
C_names = {'AF','Int','High multi','Non','Low multi'};
clusterID = zeros(height(IDRfeatures),1);
for i = 1:height(IDRfeatures)
clusterID(i,1) = find(pdist2(table2array(IDRfeatures(i,2:3)),C) == min(pdist2(table2array(IDRfeatures(i,2:3)),C)));
end
IDRfeatures(:,52) = array2table(clusterID);
IDRfeatures.Properties.VariableNames{52} = 'Cluster';
%% Test for significant differences between clusters in protein features
h_AFvsX = zeros(4,width(IDRfeatures));
p_AFvsX = zeros(4,width(IDRfeatures));
for m = 1:length(C_names)
    for j = 2:width(IDRfeatures)
        [h_AFvsX(m,j),p_AFvsX(m,j)] = ttest2(table2array(IDRfeatures(find(clusterID == 1),j)),...
            table2array(IDRfeatures(find(clusterID == m),j)),'alpha',0.05/480);  
    end
end
h_IntVsX = zeros(4,width(IDRfeatures));
p_IntVsX = zeros(4,width(IDRfeatures));
for m = 1:length(C_names)
    for j = 2:width(IDRfeatures)
        [h_IntVsX(m,j),p_IntVsX(m,j)] = ttest2(table2array(IDRfeatures(find(clusterID == 2),j)),...
            table2array(IDRfeatures(find(clusterID == m),j)),'alpha',0.05/480);
    end
end
h_HMvsX = zeros(4,width(IDRfeatures));
p_HMvsX = zeros(4,width(IDRfeatures));
for m = 1:numel(C_names)
    for j = 2:width(IDRfeatures)
        [h_HMvsX(m,j),p_HMvsX(m,j)] = ttest2(table2array(IDRfeatures(find(clusterID == 3),j)),...
            table2array(IDRfeatures(find(clusterID == m),j)),'alpha',0.05/480);
    end
end
h_NonVsX = zeros(4,width(IDRfeatures));
p_NonVsX = zeros(4,width(IDRfeatures));
for m = 1:numel(C_names)
    for j = 2:width(IDRfeatures)
        [h_NonVsX(m,j),p_NonVsX(m,j)] = ttest2(table2array(IDRfeatures(find(clusterID == 4),j)),...
            table2array(IDRfeatures(find(clusterID == m),j)),'alpha',0.05/480);
    end
end
h_LMvsX = zeros(4,width(IDRfeatures));
p_LMvsX = zeros(4,width(IDRfeatures));
for m = 1:numel(C_names)
    for j = 2:width(IDRfeatures)
        [h_LMvsX(m,j),p_LMvsX(m,j)] = ttest2(table2array(IDRfeatures(find(clusterID == 5),j)),...
            table2array(IDRfeatures(find(clusterID == m),j)),'alpha',0.05/480);
    end
end
%% Create bloxplots for features with sig differences, if plot flag is 1
if figOpt == 1
    sigCounts = sum(vertcat(h_AFvsX, h_IntVsX, h_HMvsX, h_NonVsX, h_LMvsX))>0;
    for sigNum = 1:width(IDRfeatures)
        if sigCounts(sigNum) == 1
            figure; boxplot(table2array(IDRfeatures(:,sigNum)), IDRfeatures.Cluster,'Colors','k','GroupOrder',{'4','1','2','5','3'})
            ylabel(IDRfeatures.Properties.VariableNames(sigNum) + " per TF")
            set(gca, 'xtick', 1:5); set(gca,'xticklabel',{'No effect','Active Fraction Only','Intensity Only','Low multi-tasker','High multi-tasker'})
        end
    end
end