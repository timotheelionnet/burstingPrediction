function [M, intList, g,distM]  = ...
    interactionMatrixMaker(bgData,tfList,setSelfInteractionMode)

%% build matrix
intList = unique(bgData(:));
M = zeros(numel(tfList),numel(intList));

% create hash map for TF values
tfMap = containers.Map(tfList,1:numel(tfList));
intMap = containers.Map(intList,1:numel(intList));

g = zeros(size(bgData));
for idx = 1:size(bgData,1)
    
    % fill in A x B interactions
    if isKey(tfMap,bgData(idx,1)) && isKey(intMap,bgData(idx,2))
        i = cell2mat(values(tfMap,bgData(idx,1)));
        j = cell2mat(values(intMap,bgData(idx,2)));
        M(i,j) = M(i,j) + 1;
    end
    
    % fill in B x A interactions
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