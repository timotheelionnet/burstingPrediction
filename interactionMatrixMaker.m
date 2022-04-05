function [M, intList, g,distM]  = ...
    interactionMatrixMaker(bgData,tfList,setSelfInteractionMode)

% builds an interaction matrix M from a list of interactions bgData.

%%%% INPUT
% bgData is a <number of interactions> x 2 cell array listing gene
   % symbols involved in each interaction (pairs can be redundant)
   
% tfList is the list of gene symbols defined as the union of: 
    % all TF gene symbols in tfDB, and 
    % all gene symbols for experimentally characterized genes defined in burstingData
    
% setSelfInteractionMode is a setting that decides the number of evidence assigned
    % by the script to TF self-interactions. possible values:
    % 'ones': TF self interactions are set to 1
    % 'max': TF self interactions are set to  the maximum number of
    % interactions that a TFs makes with any unique interactor


%%%%% OUTPUT
% M is the interaction matrix, where the rows are TFs and the columns their interactors 
	% the size of M is <number of elements in tfList> x <number of elements in intList>
    % the genes in rows and columns of M follow the same order as in tfList and
    % intList
    
% intList is the list of all genes involved in an interaction.

% g is a graph object: <number of interactions> x 2 array where each row
    % contains the indices i and j of the genes forming the interaction.
    % (indices correspond to the gene symbols in intList).

% distM is a square matrix with size 
    % <number of interactors in intList> x <number of interactors in intList>
    % where each cell i,j contains the shortest distance between interactor i and
    % interactor j on the interaction graph g.
    % (indices correspond to the gene symbols in intList).

%% build matrix
intList = unique(bgData(:));
M = zeros(numel(tfList),numel(intList));

% create hash map for TF values
tfMap = containers.Map(tfList,1:numel(tfList));
intMap = containers.Map(intList,1:numel(intList));

% loop through each interaction to fill the matrix 
g = zeros(size(bgData));
for idx = 1:size(bgData,1)
    
    % fill in (m,n) cell in matrix with members of the current interaction
    if isKey(tfMap,bgData(idx,1)) && isKey(intMap,bgData(idx,2))
        i = cell2mat(values(tfMap,bgData(idx,1)));
        j = cell2mat(values(intMap,bgData(idx,2)));
        M(i,j) = M(i,j) + 1;
    end
    
    % fill in (n,m) cell in matrix with members of the current interaction
    if isKey(tfMap,bgData(idx,2)) && isKey(intMap,bgData(idx,1))
        i = cell2mat(values(tfMap,bgData(idx,2)));
        j = cell2mat(values(intMap,bgData(idx,1)));
        M(i,j) = M(i,j) + 1;
        g(idx,1:2) = [i,j];
    end
    
    % fill in graph node
    if isKey(intMap,bgData(idx,1)) && isKey(intMap,bgData(idx,2))
        i = cell2mat(values(intMap,bgData(idx,1)));
        j = cell2mat(values(intMap,bgData(idx,2)));
        g(idx,1:2) = [i,j];
    end
end


%% adjust self interactions if needed
switch setSelfInteractionMode 
    case 'ones'
        for idx = 1:numel(intList)   
            % fill in A x A interactions with ones
            i = find(ismember(tfList, intList(idx)));
            if M(i,idx) == 0
                M(i,idx) = 1;
            end
            
            if sum( (g(:,1) == idx) & (g(:,2) == idx)) == 0
                g = [g;idx,idx];
            end
        end
    case 'max'
        for idx = 1:numel(intList)   
            % fill in A x A interactions with ones
            i = find(ismember(tfList, intList(idx)));
            M(i,idx) = max(M(i,:));
            
            if sum( (g(:,1) == idx) & (g(:,2) == idx)) == 0
                g = [g;idx,idx];
            end
        end
end

% remove columns of interactors that do not interact with a TF 
% this can happen if a TF X has a single interaction with a non-TF Y. The
% column corresponding to X will be all zeros since Y is not a TF there is
% no row for Y.
idx = sum(M)~=0;
M = M(:,idx);
intList = intList(idx);

%% construct graph object
g = g(prod(g,2)~=0,:);
g = graph(g(:,1),g(:,2));
distM = distances(g);

end