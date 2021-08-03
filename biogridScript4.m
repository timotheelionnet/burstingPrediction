%% prerequisites
<<<<<<< Updated upstream:biogridScript3.m
figOpt = 1;
=======
figOpt = 0;

>>>>>>> Stashed changes:biogridScript4.m
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

[M, tf,g,distM] = interactionMatrixMaker5(...
                    bgDataClean,tfList,charTFs,setSelfInteractionMode);

%% Keep only interactors that interact with both a characterized TF and a non-characterized TF
% this is needed for the training
minInteractionsPerTF = 1;
minInteractionsPerInteractor = 1;
tf2 = pruneInteractionMatrixForTraining2(M,tf,...
    minInteractionsPerTF,minInteractionsPerInteractor);

%% load Hela RNA-seq data
expressionThreshold = 0.1;
[tf3,RNAseq] = readHelaRNASeq2(RNAseqFileName,tf2,expressionThreshold);

if figOpt == 1
    [n1,x1] = hist(log(RNAseq.pTPM+1),0:.1:10);
    figure('Name',['Expression levels in Hela, threshold set to ',num2str(expressionThreshold)]);
    hold;
    plot(x1,n1/sum(n1),'DisplayName','All genes');
    [n2,x2] = hist(log(RNAseq.pTPM(ismember(RNAseq.Gene_name,intList))+1),0:.1:10);
    plot(x2,n2/sum(n2),'DisplayName','Interactors');
    xlabel('log(pTPM+1)');
    ylabel('Counts');
end

%% plot stats
if figOpt == 1
    threshMax = 200; % the following function will plot the number of 
                    % interactions/interactors as a function of the minimum
                    % number of unique interactions per interactors, for a range
                    % of 1 to threshMax
    threshInteractions = dispInteractions2(M,burstingData, tf3,threshMax,...
    geneSymbolToDisplay1,geneSymbolToDisplay2);
end
%% normalize interaction matrix
normByRows = 1;
normByCols = 1;
M3norm = normalizeInteractionsMatrix(M3,normByRows,normByCols);

%% cluster interaction matrix
dendrogramInteractionThreshold = 50;
if figOpt == 1
<<<<<<< Updated upstream:biogridScript3.m
cg1 = clustergram(M3norm(sum(M3'>0)>10,:),'Colormap',redbluecmap,'DisplayRange',10,...
    'RowPDist','correlation','Symmetric',false,'Rowlabels',tfList2(sum(M3'>0)>10),...
    'ColumnLabels',intList3);
%%
cg2 = clustergram(M3norm,'Colormap',redbluecmap,'DisplayRange',10,...
    'Symmetric',false,'Rowlabels',tfList2,...
    'ColumnLabels',intList3);
c1 = clusterGroup(cg2,501,1);
=======
    cg1 = clustergram(M3norm(sum(M3'>0)>dendrogramInteractionThreshold,:),'Colormap',redbluecmap,'DisplayRange',10,...
        'RowPDist','correlation','Symmetric',false,'Rowlabels',tfList2(sum(M3'>0)>dendrogramInteractionThreshold),...
        'ColumnLabels',intList3);
    
%     cg2 = clustergram(M3,'Colormap',redbluecmap,'DisplayRange',10,...
%         'Symmetric',false,'Rowlabels',tfList2,...
%         'ColumnLabels',intList3);
%     c1 = clusterGroup(cg2,501,1);
>>>>>>> Stashed changes:biogridScript4.m
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
