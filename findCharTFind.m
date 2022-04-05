function [charTFind, burstingData2] = findCharTFind(charTFs, tfList, burstingData)
% Prune experimental data ensuring that that TFS 
% that are experimentally characterized factors (i.e. rows of burstingData)
% are also present in the list of all TFs (tfList2). 
% charTFind are the indices of the retained factors relative to the gene symbols in tfList 

%%%%% INPUT
% charTFs: list of all gene symbols for experimentally characterized genes defined in burstingData
% tfList: list of gene symbols for the rows of the interaction matrix (TFs)
% burstingData: MATLAB table with columns TFnames, activity, intensity (CRISPRburst results)  
        % the only column used here is 'TFnames'

%%%%% OUTPUT
% charTFind: indices of experimentally characterized genes (i.e. present in
        % burstingData), relative to the list of gene symbols in tfList.
% burstingData2: burstingData where the rows corresponding to genes absent
    % from tfList have been removed.

%%

% for each experimentally characterized factor (from charTFs list), 
% collect the index of its location in tfList (if present)
charTFind = [];
for v = 1:length(charTFs)
    charTFind = [charTFind, strmatch(charTFs(v), tfList, 'exact')];
end

% find out which experimentally characterized TFs are also present in the list
% of TFs, and collect their indices
retained = [];
for f = 1:numel(charTFind)
    retained = [retained, strmatch(tfList(charTFind(f)), burstingData.TFnames,'exact')];
end

% prune experimental data to keep only the characterized genes that have a match in the TF list.
burstingData2 = burstingData(retained,:);

end