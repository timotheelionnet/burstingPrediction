%% prerequisites
figOpt = 0;
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
minInteractionsPerInteractor = 5;
[M2,intList2,isInteractorInBothLists,tfList2,lonelyTFs,~] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData,...
    minInteractionsPerTF,minInteractionsPerInteractor);

distM2 = distM(isInteractorInBothLists,isInteractorInBothLists);

%% load Hela RNA-seq data
expressionThreshold = 0.1;
[M3,intList3,RNAseq] = readHelaRNASeq(RNAseqFileName,intList2,expressionThreshold, M2);

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
if figOpt == 1
cg1 = clustergram(M3(sum(M2'>0)>30,:),'Colormap',redbluecmap,'DisplayRange',10,...
    'RowPDist','correlation','Symmetric',false,'Rowlabels',tfList2(sum(M2'>0)>30),...
    'ColumnLabels',intList3);
%%
cg2 = clustergram(M3,'Colormap',redbluecmap,'DisplayRange',10,...
    'Symmetric',false,'Rowlabels',tfList2,...
    'ColumnLabels',intList3);
c1 = clusterGroup(cg2,501,1);
end
%% Find indices of characterized TFs in list of all TFs and subset burstingData for training
[charTFind, burstingData2] = findCharTFind(charTFs, tfList2, burstingData);

%% Train and predict
%Threshold setting is the minimum number an interactor must make to be included in the prediction
%Insert interaction threshold for Activity here
A = 30;
%Insert interaction threshold for Intensity here
I = 15;
[predScores,betaTrainAct,betaTrainInt] = updatingModelwip ...
    (M3,M3norm,tfList2,charTFind,burstingData2,A,I);

%% Plot predictions
plotPredictions(M3,predScores,burstingData2,tfList2,charTFind,A,I)
%% 
toc   
