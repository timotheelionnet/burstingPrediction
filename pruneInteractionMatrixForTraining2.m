function tf2 = pruneInteractionMatrixForTraining2(M,tf,...
    minInteractionsPerTF,minInteractionsPerInteractor)

% remove interactors that do not make interactions with the characterized
% and non-characterized group

%% check that each interactor makes at least one interaction with a
% characterized TF and at least one interaction with a non-characterized
% TF.
isChar = cell2mat(values(tf.map,tf.char));
isNonChar = cell2mat(values(tf.map,tf.nonChar));

isInteractorInBothLists = sum(M(isChar,:))>0 & sum(M(isNonChar,:))>0;

tf2 = tf;
tf2.trainInt = tf.int(isInteractorInBothLists);

%% remove interactors that have fewer than minInteractionsPerInteractor
hasInteractorEnoughInteractions = sum(M) >= minInteractionsPerInteractor;
tf2.trainInt = intersect(tf2.trainInt, tf.int(hasInteractorEnoughInteractions));

%% remove TFs that do not have any interactors (or fewer interactors than minInteractionsPerTF)
isTF = cell2mat(values(tf.map,tf.allTFs));
hasTFEnoughInteractions = sum(M(isTF,:),2) >= minInteractionsPerTF;
tf2.trainTF = tf.allTFs(hasTFEnoughInteractions);


end