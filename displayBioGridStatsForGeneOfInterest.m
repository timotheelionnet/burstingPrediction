function displayBioGridStatsForGeneOfInterest(geneSymbolToDisplay,bgData)

%% 
ui1 = unique(bgData(ismember(bgData(:,1),geneSymbolToDisplay),2)) ;
ui2 = unique(bgData(ismember(bgData(:,2),geneSymbolToDisplay),1)) ;
ui = numel(unique(union(ui1,ui2)));
n1 = numel( bgData(ismember(bgData(:,1),geneSymbolToDisplay),2) );
n2 = numel( bgData(ismember(bgData(:,2),geneSymbolToDisplay),1) );
disp([geneSymbolToDisplay,': ',num2str(ui),' total unique interactors; ',num2str(n1+n2),...
    ' total interactions. Interactors/interactions as A: ',num2str(numel(ui1)),...
    '/',num2str(n1),'; Interactors/interactions as B: ',num2str(numel(ui2)),'/',num2str(n2)]);
% b = sortrows(bgData(ismember(bgData(:,1),geneSymbolToDisplay) ...
%     | ismember(bgData(:,2),geneSymbolToDisplay),:),2);
% b = sortrows(b,1);