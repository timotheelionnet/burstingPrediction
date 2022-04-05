function [bgData,tfList,charTFs] = pruneInteractionList(bgData,tfDB,burstingData)
% From a list of genome-wide interaction pairs (bgData), 
% this function keeps only interactions where at least one partner is a TF

%%%%% INPUT: 
% bgData is a <number of interactions> x 2 cell array listing gene
            % symbols involved in each interaction (pairs can be redundant)
% tfDB is a MATLAB table containing TFs as rows and their features as columns            
        % the only column used here is 'HGNCSymbol'
% burstingData is a MATLAB table with columns TFnames, activity, intensity (CRISPRburst results)  
        % the only column used here is 'TFnames'
        
%%%%% OUTPUT:
% bgData is the pruned list of interactions
% tfList is the list of gene symbols defined as the union of: 
    % all TF gene symbols in tfDB, and 
    % all gene symbols for experimentally characterized genes defined in burstingData
% charTFs is the list of all gene symbols for experimentally characterized genes defined in burstingData
    
%%

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