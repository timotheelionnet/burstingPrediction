function plotPredictions(predScores,burstingData,tfList,charTFind,exampleTFInd,characterizedInd)

% plot kinetic predictions as a scatter plot

%%%%% INPUT
% predScores: predicted kinetics for all TFs ; array with size <number of TFs> x 3 where:
        % col 1 is the predicted active fraction for each TF, 
        % col 2 is the predicted intensity for each TF
        % col 3 is the TF index corresponding to the gene symbol list tfList
        % row indices should correspond to gene symbols in tfList
        
% burstingData: MATLAB table holding measured kinetic data 
        % with columns TFnames, activity, intensity (CRISPRburst results)  
        
% tfList: list of gene symbols for the rows of the kinetics in predScores

% charTFind: indices mapping experimentally characterized genes (i.e. rows of
        % burstingData) to the list of gene symbols in tfList.
        
% exampleTFInd: indices (relative to tfList) for arbitrary TF(s) chosen as examples

% characterizedInd: indices of the characterized TFs to plot (relative to rows of burstingData2 matrix)

%%
% subset of the TFs that have not been characterized
predictIndex = setdiff(1:numel(tfList),charTFind);

figure('Name','TF kinetic predictions');

% plot predicted kinetics for the all non-experimentally characterized TFs
% (semi transparent gray)
scatter(predScores(predictIndex,1),predScores(predictIndex,2), 'filled', ...
    'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', '#262521', 'MarkerEdgeAlpha',0.0, ...
    'SizeData', 20, 'MarkerFaceColor', '#262521'); 
hold on;

% plot experimentally measured kinetics for the measured TFs
% (~solid orange)
scatter(burstingData.activity, burstingData.intensity,'filled',...
    'MarkerFaceAlpha', 0.8,'MarkerEdgeAlpha',0.8,'MarkerEdgeColor','#fc9403', ...
    'SizeData', 20, 'MarkerFaceColor', '#fc9403'); 

% label the data points for the measured TFs with their gene symbols
labelpoints(burstingData.activity(characterizedInd),...
    burstingData.intensity(characterizedInd), ...
    tfList(charTFind(characterizedInd)), 'N', 0.01); 

% plot the predicted scores for the example TFs chosen
% (~solid blue)
scatter(predScores(exampleTFInd,1),predScores(exampleTFInd,2), 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', '#2403fc', 'MarkerEdgeAlpha',0.8, ...
    'SizeData', 20, 'MarkerFaceColor', '#2403fc'); 

% label the data points for the example TFs with their gene symbols
labelpoints(predScores(exampleTFInd,1),...
    predScores(exampleTFInd,2), ...
    tfList(exampleTFInd), 'N', 0.01); hold on;

% plot dashed lines x = 1 and y=1 as reference
plot([0,5],[1,1], 'LineStyle', ':', 'Color', 'k'); 
plot([1,1],[0,2], 'LineStyle', ':', 'Color', 'k');

% set plot size, aspect ratio, axis legends etc
xlabel('Predicted Active Fraction Effect');
ylabel('Predicted Intensity Effect');
xlim([0.0 3.5]); 
ylim([0.5 1.8]); 
box off; grid off; 
set(gca,'fontsize', 12);
pbaspect([1.75 1 1]);

end