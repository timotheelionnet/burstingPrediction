%% prerequisites
figOpt = 1;
% 1) BIOGRID database file: we used the .tab3 formatted file from 
% https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.199/
% set name of the downladed file here:
bioGridFileName = 'BIOGRID-ALL-4.4.199.tab3.txt';

% 2) Bursting data:
burstingDataFileName = 'burstingData.txt';

% 3) Full human TF database, from Lambert et al (PMID:29425488) 
% downloaded as txt format from http://humantfs.ccbr.utoronto.ca/download.php
humanTFDataBaseFileName = 'DatabaseExtract_v_1.01.txt';

% 4) RNA-seq data for HeLa cell line, from
RNAseqFileName = "hela_RNAseq.xlsx";

tic
%% load the list of interactions between all TFs
geneSymbolToDisplay1 = 'EP300';   % a first arbitrary gene for which the stats 
                                % will be displayed in the command line.
geneSymbolToDisplay2 = 'BRD4';   % a second arbitrary gene for which the stats 
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
    pruneInteractionList3(bgData,tfDB,burstingData);

%% create interaction matrix
setSelfInteractionMode = 'max'; % this setting sets the number of evidence of a TF self interaction
                                    % to the maximum number of evidence
                                    % that TF has for any factor

[M, intList,g,distM] = interactionMatrixMaker4(...
                    bgDataClean,tfList,setSelfInteractionMode);

%% Keep only interactors that interact with both a characterized TF and a non-characterized TF
% this is needed for the training
minInteractionsPerTF = 1;
minInteractionsPerInteractor = 1;
[M2,intList2,isInteractorInBothLists,tfList2,lonelyTFs,~] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData,...
    minInteractionsPerTF,minInteractionsPerInteractor);

distM2 = distM(isInteractorInBothLists,isInteractorInBothLists);

%% load Hela RNA-seq data
expressionThreshold = 0.1;
[M3,intList3,RNAseq] = readHelaRNASeq(RNAseqFileName,intList2,expressionThreshold, M2);

map2 = containers.Map(intList2,1:numel(intList2));
idx3 = cell2mat(values(map2,intList3));
distM3 = distM2(idx3,idx3);

%% plot stats
if figOpt == 1
    threshMax = 200; % the following function will plot the number of 
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
minInteractionsPerTFClusterGram = 30;
if figOpt == 1
    
    % dendrogram of all interactors by distance - gets pretty large
    idxTF3 = sum(M3'>0)>minInteractionsPerTFClusterGram;
    map3 = containers.Map(intList3,1:numel(intList3));
    int3 = setdiff(intList3, tfList2(~idxTF3));
    idxM3 = cell2mat(values(map3,int3));
%     invDist3 = 1./distM3;
%     invDist3(isinf(invDist3)) = 2;
%     invDist3 = log(invDist3+1);
    cg1 = clustergram(-distM3(idxM3,idxM3),...
        'Colormap',redbluecmap,'DisplayRange',0.5,...
        'RowPDist','correlation',...
        'ColumnPDist','correlation','Symmetric',false,...
        'Rowlabels',int3,...
        'ColumnLabels',int3);
    addTitle(cg1,'all interactors clustered by interaction distance');
    
    % same as previous but use raw interactions (0/1) rather than distance
    x = zscore(zscore(distM3(idxM3,idxM3) <= 1')');
    x = (x+x')/2; % enforcing that the matrix stays symmetric
    cg1b = clustergram(x,...
        'Colormap',redbluecmap,'DisplayRange',1,...
        'RowPDist','correlation',...
        'ColumnPDist','correlation','Symmetric',false,...
        'Rowlabels',int3,...
        'ColumnLabels',int3);
    addTitle(cg1b,'all interactors clustered by raw interactions');
    
    % restricting the dendrogram to TFs only 
    idxTF4 = find(ismember(intList3,tfList2));
%     invDist3 = 1./distM3;
%     invDist3(isinf(invDist3)) = 2;
%     invDist3 = log(invDist3+1);
    
    cg2 = clustergram(-distM3(idxTF4,idxTF4),...
         'Colormap',redbluecmap,...
         'DisplayRange',4,...
         'RowPDist','correlation',...
         'ColumnPDist','correlation','Symmetric',false,...
         'Rowlabels',intList3(idxTF4),...
         'linkage','complete',...
         'ColumnLabels',intList3(idxTF4));
     
    addTitle(cg2,'TFs only clustered by interaction distance');
    
    % same as previous but use raw interactions (0/1) rather than distance
    x = zscore(zscore(distM3(idxTF4,idxTF4) <= 1')');
    x = (x+x')/2; % enforcing that the matrix stays symmetric
    cg2b = clustergram(x,...
         'Colormap',redbluecmap,...
         'DisplayRange',4,...
         'RowPDist','correlation',...
         'ColumnPDist','correlation','Symmetric',false,...
         'Rowlabels',intList3(idxTF4),...
         'linkage','complete',...
         'ColumnLabels',intList3(idxTF4));
    
    addTitle(cg2b,'TFs only clustered by raw interactions');

end


%% Find indices of characterized TFs in list of all TFs and subset burstingData for training
[charTFind, burstingData2] = findCharTFind(charTFs, tfList2, burstingData);

%% Train and predict
%Threshold setting is the minimum number an interactor must make to be included in the prediction
%Insert interaction threshold for Activity here
A = 60;
%Insert interaction threshold for Intensity here
I = 25;
[predScores,betaTrainAct,betaTrainInt] = updatingModelwip ...
    (M3,M3norm,tfList2,charTFind,burstingData2,A,I);

%% Plot predictions
plotPredictions(M3,predScores,burstingData2,tfList2,charTFind,A,I)

%% Find top predictors using beta values calculated by model 
[topPredictors_activity, topPredictors_intensity] ...
    = findPredictors(betaTrainAct, betaTrainInt, M3, intList3, A, I);

%%
toc   
