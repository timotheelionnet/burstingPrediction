function [M, tf, g,distM]  = ...
    interactionMatrixMaker5(bgData,tfList,charTFs,setSelfInteractionMode)

% create hash map for TF values
intList = unique(bgData(:));
fullList = union(tfList,union(charTFs,intList));

tf = struct();
tf.map = containers.Map(fullList,1:numel(fullList));
tf.allTFs = tfList;
tf.char = charTFs;
tf.nonChar = setdiff(tfList,charTFs);
tf.int = fullList;

%% build matrix
M = zeros(numel(tf.int));
g = zeros(size(bgData));
for idx = 1:size(bgData,1)
    
    % fill in interactions
    if isKey(tf.map,bgData(idx,1)) && isKey(tf.map,bgData(idx,2))
        i = cell2mat(values(tf.map,bgData(idx,1)));
        j = cell2mat(values(tf.map,bgData(idx,2)));
        M(i,j) = M(i,j) + 1;
        M(j,i) = M(j,i) + 1;
        g(idx,1:2) = [i,j];
    end
end
g = g(prod(g,2)~=0,:);
g = graph(g(:,1),g(:,2));
distM = distances(g);

%% adjust self interactions if needed
switch setSelfInteractionMode 
    case 'ones'
        for idx = 1:numel(numel(tf.int))   
            % fill in A x A interactions with ones
            if M(idx,idx) == 0
                M(idx,idx) = 1;
            end
        end
    case 'max'
        for idx = 1:numel(intList)   
            % fill in A x A interactions with maximum interactions that
            % factor makes
            M(idx,idx) = max(M(idx,:));
        end
end

end