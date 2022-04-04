function plotPredictions(M,predScores,burstingData,tfList,charTFind,A,I,Ind)

InteractionIndexStoreAct = find(sum(M)>=A);
InteractionIndexStoreInt = find(sum(M)>=I);
predictIndex = setdiff(1:numel(tfList),charTFind);
figure;
oi = [1:72];
scatter(predScores(predictIndex,1),predScores(predictIndex,2), 'filled', 'MarkerFaceAlpha', 0.1, ...
    'MarkerEdgeColor', '#262521', 'MarkerEdgeAlpha',0.0, 'SizeData', 20, 'MarkerFaceColor', '#262521'); hold on;

scatter(burstingData.activity, burstingData.intensity,'filled','MarkerFaceAlpha', 0.8,'MarkerEdgeAlpha',0.8,'MarkerEdgeColor','#fc9403', ...
'SizeData', 20, 'MarkerFaceColor', '#fc9403'); hold on;
labelpoints(burstingData.activity(oi),burstingData.intensity(oi), tfList(charTFind(oi)), 'N', 0.01); hold on;

scatter(predScores(Ind,1),predScores(Ind,2), 'filled', 'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeColor', '#2403fc', 'MarkerEdgeAlpha',0.8, 'SizeData', 20, 'MarkerFaceColor', '#2403fc'); hold on;
labelpoints(predScores(Ind,1),predScores(Ind,2), tfList(Ind), 'N', 0.01); hold on;

xlabel(['Predicted Active Fraction Effect']);
ylabel(['Predicted Intensity Effect']);
plot([0,5],[1,1], 'LineStyle', ':', 'Color', 'k'); plot([1,1],[0,2], 'LineStyle', ':', 'Color', 'k');
xlim([0.0 3.5]); ylim([0.5 1.8]); box off; grid off; set(gca,'fontsize', 12);
pbaspect([1.75 1 1])

end