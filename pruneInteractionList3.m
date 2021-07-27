function [bgData,tfList] = pruneInteractionList3(bgData,tfDB,burstingData)

%% keep only interactions where at least one partner is a TF

% list of names of characterized TFs
charTFs = unique(burstingData.TFnames(:));

% list of names of all TFs from human TF database
nonCharTFs = unique(tfDB.HGNCSymbol(:));

% combine all TFs into a single list
tfList = union(charTFs, nonCharTFs);

% keep only interactions that involve at least one TF
bgData = bgData(...
    ismember(bgData(:,1),tfList) | ismember(bgData(:,2),tfList) ,: ); 

% sort rows
bgData = sortrows(bgData,[1,2]);

disp(['There are ',num2str(numel(union(bgData(:,1),bgData(:,2)))),...
    ' interactors that interact with a TF.']);

end