%% prerequisites

% 1) BIOGRID database file: we used the .tab3 formatted file from 
% https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.199/
% set name of the downladed file here:
bioGridFileName = 'BIOGRID-ALL-4.4.199.tab3.txt';

% 2) Bursting data:
burstingDataFileName = 'burstingData.txt';

% 3) Full human TF database, from Lambert et al (PMID:29425488) 
% downloaded as txt format from http://humantfs.ccbr.utoronto.ca/download.php
humanTFDataBaseFileName = 'DatabaseExtract_v_1.01.txt';

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
[bgDataClean,tfList] = ...
    pruneInteractionList3(bgData,tfDB,burstingData);

%% create interaction matrix
setSelfInteractionMode = 'max'; % this setting sets the number of evidence of a TF self interaction
                                    % to the maximum number of evidence
                                    % that TF has for any factor

[M, intList,g,distM] = interactionMatrixMaker4(...
                    bgDataClean,tfList,setSelfInteractionMode);

%% Keep only interactors that interact with both a characterized TF and a non-characterized TF
% this is needed for the training
[M2,intList2,isInteractorInBothLists] = ...
    pruneInteractionMatrixForTraining(M,intList,tfList,burstingData); 

distM2 = distM(isInteractorInBothLists,isInteractorInBothLists);
%% plot stats
threshMax = 200; % the following function will plot the number of 
                    % interactions/interactors as a function of the minimum
                    % number of unique interactions per interactors, for a range
                    % of 1 to threshMax
threshInteractions = dispInteractions(M2,burstingData, tfList,intList2,threshMax,...
    geneSymbolToDisplay1,geneSymbolToDisplay2);

%% truncate interactor list by setting a threshold number of interactions
minInteractions = 30;
[M2,intList2] = cutoffAndNormalizeInteractionsMatrix(M,intList,minInteractions);

%% cluster interaction matrix

cg = clustergram(M2(sum(M'>0)>30,:),'Colormap',redbluecmap,'DisplayRange',10,...
    'RowPDist','correlation','Symmetric',false,'Rowlabels',tfList(sum(M'>0)>30),...
    'ColumnLabels',intList2);
%%
cg = clustergram(M2,'Colormap',redbluecmap,'DisplayRange',10,...
    'Symmetric',false,'Rowlabels',tfList,...
    'ColumnLabels',intList2);
%%
% cg = clustergram(M,'Colormap',redbluecmap,'DisplayRange',10,...
%     'Symmetric',false,'Rowlabels',tfList,...
%     'ColumnLabels',intList);
%%
%c1 = clusterGroup(cg,501,1);
%%
%charTF_IntList_from_allTF_IntList;

%% 
%charTF_train_allTF_predict;

%% 
toc   