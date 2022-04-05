%% prerequisites
figOpt = 1;
% 1) BIOGRID database file: we used the .tab3 formatted file from
% https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.207/
% set name of the downloaded file here:
bioGridFileName = '../burstingPredictionData/BIOGRID-ALL-4.4.207.tab3.txt';

% 2) Bursting data:
burstingDataFileName = 'burstingData.txt';

% 3) Full human TF database, from Lambert et al (PMID:29425488)
% downloaded as txt format from http://humantfs.ccbr.utoronto.ca/download.php
humanTFDataBaseFileName = 'DatabaseExtract_v_1.01.txt';

% 4) RNA-seq data for HeLa cell line, from
RNAseqFileName = "hela_RNAseq.xlsx";

addpath('labelpoints.m');

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

%% keep only interactions that involved at least one TF
[bgDataClean,tfList,charTFs] = ...
    pruneInteractionList(bgData,tfDB,burstingData);

%% create interaction matrix from list of interactions
setSelfInteractionMode = 'max'; 
% this setting sets the number of evidence of a TF self-interaction
% to the maximum number of evidence that TF has for any factor

% builds the interaction matrix M from the list of interactions bgDataClean.
[M, intList,g,distM] = interactionMatrixMaker(...
    bgDataClean,tfList,setSelfInteractionMode);

%% Keep only interactors that interact with both a characterized TF and a non-characterized TF
% this is needed for the training
minInteractionsPerTF = 5;   % each TF must have at least 5 interactions 
                            % to be used in the training set
minInteractionsPerInteractor = 5; % each Interactor must have at least 5 
                              % interactions to be used in the training set
                              
% M2 is the pruned interaction matrix, with rows corresponding to tfList2 and columns to intList2                              
[M2,intList2,isInteractorInBothLists,tfList2,lonelyTFs,~] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData,...
    minInteractionsPerTF,minInteractionsPerInteractor);

% distM2 is a square matrix that contains the distance between gene i and gene j along the
% interaction network, rows/columns indices correspond to distGenes2.
distGenes2 = intList(isInteractorInBothLists);
distM2 = distM(isInteractorInBothLists,isInteractorInBothLists);

%% load Hela RNA-seq data in order to select only interactors expressed in HeLa
expressionThreshold = 0.1; % minimum TPM value in HeLa RNA-seq data for 
                           % a gene to be considered as expressed in HeLa

% prune interaction matrix M2 to exclude interactors not expressed in HeLa.                           
[M3,intList3,RNAseq] = readHelaRNASeq(...
    RNAseqFileName,intList2,expressionThreshold, M2, figOpt);

% distM3 is a square matrix that contains the distance between gene i and gene j along the
% interaction network, rows/columns indices correspond to intList3.
map2 = containers.Map(intList2,1:numel(intList2));
idx3 = cell2mat(values(map2,intList3));
distM3 = distM2(idx3,idx3);

%% plot statistics for the number of interactions (optional)
if figOpt == 1
    threshMax = 500; 
    
    % plot the number of interactions/interactors as a function of the 
    %minimum number of unique interactions per interactors, for a range of 1 to threshMax
    threshInteractions = dispInteractions(M3,burstingData, tfList2,intList3,threshMax,...
        geneSymbolToDisplay1,geneSymbolToDisplay2);
end

%% normalize interaction matrix
normByRows = 1;
normByCols = 1;
M3norm = normalizeInteractionsMatrix(M3,normByRows,normByCols);

%% cluster interaction matrix and distance matrix and display as clustergrams (optional)
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

%% Prune experimental data for training
% ensuring that that TFS that are experimentally characterized (rows of burstingData)
% are also present in the list of all TFs (tfList2). 
% charTFind are the indices of the retained factors (rows of burstingData2),
% relative to the gene symbols in tfList2 
[charTFind, burstingData2] = findCharTFind(charTFs, tfList2, burstingData);

%% Train the regression model that predicts TF kinetics based on interaction data
% for a range of interaction thresholds

% interactors to use in the regression are retained if they make more than
% a certain threshold number of total interactions. The following code section 
% loops through a range of those threshold interaction numbers.
range = 10; % the number of different threshold values to loop through.
Amin = 235; % the starting threshold value for the active fraction model
Imin = 25; % the starting threshold value for the intensity model

nCompA = 4; % number of partial least square components to use in regression model for active fraction
nCompI = 6; % number of partial least square components to use in regression model for intensity

RollingScores = {};
for th = 1:range
    minInteractionsA(th) = Amin+th; % minimum number of interactions for an interactor 
                                % to be used in the active fraction model
                                
    minInteractionsI(th) = Imin+th; % minimum number of interactions for an interactor 
                                % to be used in the intensity model
    
    % train model, output predictive weights for retained interactors (beta) 
    % as well as predicted kinetic scores for all TFs in M (predScores)                             
    [predScores,betaTrainAct{th},betaTrainInt{th}] ...
        = TFregress(M3,M3norm,tfList2,charTFind,burstingData2,...
        minInteractionsA(th),minInteractionsI(th),nCompA,nCompI);

    RollingScores(th,1:2) = {minInteractionsI(th),predScores};
end

%% Plot kinetic predictions of the model for all TFs

% indices (relative to tfList2) for the two TFs chosen as examples:
exampleTFInd = [ strmatch(geneSymbolToDisplay1, tfList2,'exact'), ...
    strmatch(geneSymbolToDisplay2, tfList2,'exact')];

% indices of the characterized TFs to plot (relative to rows of burstingData2 matrix):
characterizedInd = 1:72; 

threshIndex = 5; % choosing the threshold setting for the predictions to plot

% plot predictions
plotPredictions(RollingScores{threshIndex,2},burstingData2,tfList2,...
    charTFind,exampleTFInd,characterizedInd)

%% Find interactor genes that are top predictors using beta values calculated by model

[topPredictorsActivity, topPredictorsIntensity] ...
    = findPredictors(...
    betaTrainAct{threshIndex}, betaTrainInt{threshIndex},...
    M3, intList3, ...
    minInteractionsA(threshIndex), minInteractionsI(threshIndex));

%%
toc
