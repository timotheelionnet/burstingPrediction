function plotPredictions(M,predScores,burstingData,tfList,charTFind,A,I)
InteractionIndexStoreAct = find(sum(M)>=A);
InteractionIndexStoreInt = find(sum(M)>=I);
figure;

scatter(predScores(:,1),predScores(:,2), 'filled', 'MarkerFaceAlpha', 0.2, ...
'MarkerEdgeColor', '#262521', 'MarkerEdgeAlpha',0.2, 'SizeData', 20, 'MarkerFaceColor', '#262521'); hold on;

%scatter(predAllScores(libTFs_ind,1),predAllScores(libTFs_ind,2), 'filled', 'MarkerFaceAlpha', 0.6, ...
%'MarkerEdgeColor', '#2b5cab', 'MarkerEdgeAlpha',0.6, 'SizeData', 20, 'MarkerFaceColor', '#2b5cab'); hold on;
%labelpoints(predAllScores(libTFs_ind,1),predAllScores(libTFs_ind,2),tfList2(predAllScores(libTFs_ind,3)),'N',0.01); hold on;

scatter(burstingData.activity, burstingData.intensity,'filled','MarkerFaceAlpha', 0.6,'MarkerEdgeAlpha',0.8,'MarkerEdgeColor','#fc9403', ...
'SizeData', 20, 'MarkerFaceColor', '#fc9403'); hold on;
labelpoints(burstingData.activity,burstingData.intensity, tfList(charTFind), 'N', 0.01); hold on;

predScores = sortrows(predScores,1,'descend');
labelpoints(predScores(1:20,1),predScores(1:20,2),tfList(predScores(1:20,3)),'N',0.01); hold on;
predScores = sortrows(predScores,2,'descend');
labelpoints(predScores(1:20,1),predScores(1:20,2),tfList(predScores(1:20,3)),'N',0.01); hold on;

xlabel(['Predicted Activity Effect with ' num2str(numel(InteractionIndexStoreAct)) ' interactors']);
ylabel(['Predicted Intensity Effect with ' num2str(numel(InteractionIndexStoreInt)) ' interactors']);
plot([0,3.5],[1,1], 'LineStyle', ':', 'Color', 'k'); plot([1,1],[0,2], 'LineStyle', ':', 'Color', 'k');
xlim([0 inf]); ylim([0 inf]); box on; grid on; set(gca,'fontsize', 12);

end