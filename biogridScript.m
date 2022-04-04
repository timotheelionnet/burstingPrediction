%% prerequisites
figOpt = 0;
% 1) BIOGRID database file: we used the .tab3 formatted file from
% https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.207/
% set name of the downloaded file here:
bioGridFileName = 'BIOGRID-ALL-4.4.207.tab3.txt';

% 2) Bursting data:
burstingDataFileName = 'burstingData.txt';

% 3) Full human TF database, from Lambert et al (PMID:29425488)
% downloaded as txt format from http://humantfs.ccbr.utoronto.ca/download.php
humanTFDataBaseFileName = 'DatabaseExtract_v_1.01.txt';

% 4) RNA-seq data for HeLa cell line, from
RNAseqFileName = "hela_RNAseq.xlsx";

tic
%% load the list of interactions between all TFs
geneSymbolToDisplay1 = 'KLF4';   % a first arbitrary gene for which the stats
% will be displayed in the command line.
geneSymbolToDisplay2 = 'GATA1';   % a second arbitrary gene for which the stats
% will be displayed in the command line.
bgData = readBioGridFile(bioGridFileName,geneSymbolToDisplay1);

displayBioGridStatsForGeneOfInterest(geneSymbolToDisplay1,bgData);
displayBioGridStatsForGeneOfInterest(geneSymbolToDisplay2,bgData);
disp(' ');

%% load bursting data
burstingData = readtable(burstingDataFileName);
disp(['Loaded bursting data from ',...
    num2str(size(burstingData,1)) ,' genes.'] );

%% load TF database
tfDB = readtable(humanTFDataBaseFileName);

%% keep only interactors that interact with one of the TFs
[bgDataClean,tfList,charTFs] = ...
    pruneInteractionList(bgData,tfDB,burstingData);
%tat opt here
%% create interaction matrix
setSelfInteractionMode = 'max'; % this setting sets the number of evidence of a TF self interaction
% to the maximum number of evidence
% that TF has for any factor

[M, intList,g,distM] = interactionMatrixMaker(...
    bgDataClean,tfList,setSelfInteractionMode);

%% Keep only interactors that interact with both a characterized TF and a non-characterized TF
% this is needed for the training
minInteractionsPerTF = 5;
minInteractionsPerInteractor = 5;
[M2,intList2,isInteractorInBothLists,tfList2,lonelyTFs,~] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData,...
    minInteractionsPerTF,minInteractionsPerInteractor);

distM2 = distM(isInteractorInBothLists,isInteractorInBothLists);

%% load Hela RNA-seq data
expressionThreshold = 0.1;
[M3,intList3,RNAseq] = readHelaRNASeq(RNAseqFileName,intList2,expressionThreshold, M2, figOpt);

map2 = containers.Map(intList2,1:numel(intList2));
idx3 = cell2mat(values(map2,intList3));
distM3 = distM2(idx3,idx3);

%% plot stats
if figOpt == 1
    threshMax = 500; % the following function will plot the number of
    % interactions/interactors as a function of the minimum
    % number of unique interactions per interactors, for a range
    % of 1 to threshMax
    threshInteractions = dispInteractions(M3,burstingData, tfList2,intList3,threshMax,...
        geneSymbolToDisplay1,geneSymbolToDisplay2);
end
%% normalize interaction matrix
normByRows = 1;
normByCols = 1;
M3norm = normalizeInteractionsMatrix(M3,normByRows,normByCols);

%% cluster interaction matrix
minInteractionsPerTFClusterGram = 10;
if figOpt == 1
    idxTF3 = sum(M3'>0)>minInteractionsPerTFClusterGram;
        cg1 = clustergram(M3norm(idxTF3,:),...
            'Colormap',redbluecmap,'DisplayRange',10,...
            'RowPDist','correlation','Symmetric',false,...
            'Rowlabels',tfList2(idxTF3),...
            'ColumnLabels',intList3);


    map3 = containers.Map(intList3,1:numel(intList3));
    int3 = setdiff(intList3, tfList2(~idxTF3));
    idxM3 = cell2mat(values(map3,int3));
    invDist3 = 1./distM3;
    invDist3(isinf(invDist3)) = 2;
    invDist3 = log(invDist3+1);
    cg2 = clustergram(invDist3(idxM3,idxM3),...
        'Colormap',redbluecmap,'DisplayRange',0.6,...
        'RowPDist','euclidean','Symmetric',false,...
        'Rowlabels',int3,...
        'ColumnLabels',int3, 'LogTrans', true);

end
%% Find indices of characterized TFs in list of all TFs and subset burstingData for training
[charTFind, burstingData2] = findCharTFind(charTFs, tfList2, burstingData);

%% Train and predict
%Threshold setting is the minimum number an interactor must make to be included in the prediction
RollingScores = {};
range = 10;
for th = 1:range
    %Insert interaction threshold for Activity here
    A = 235+th;
    %Insert interaction threshold for Intensity here
    I = 25+th;
    [predScores,betaTrainAct,betaTrainInt] ...
    = TFregress(M3,M3norm,tfList2,charTFind,burstingData2,A,I);

    RollingScores(th,1:2) = {I,predScores};
end
%% Plot predictions
% TF of interest example
ind = [strmatch(geneSymbolToDisplay1, tfList2,'exact'),strmatch(geneSymbolToDisplay2, tfList2,'exact')];
plotPredictions(M3,RollingScores{5,2},burstingData2,tfList2,charTFind,A,I,ind)

%% Find top predictors using beta values calculated by model
[topPredictors_activity, topPredictors_intensity] ...
    = findPredictors(betaTrainAct, betaTrainInt, M3, intList3, A, I);

%%
toc
