function [M2,intList2,RNAseq] = readHelaRNASeq(RNAseqFileName,intList,expressionThreshold, M, figOpt)
% prune interaction matrix M to exclude interactors (columns) that are 
% not expressed in HeLa based on RNAseq data and an expression threshold.  

%%%%% INPUT
% RNAseqFileName: path to the excel spreadsheet holding the RNA-seq TPM data
% intList: list of interactors gene symbols corresponding to the columns of M.
% expressionThreshold: minimum TPM value for a gene to be considered
        % expressed.
% M: interaction matrix, where the columns are the interactors in intList.        
% figOpt: set this flag to 1 to display distribution of TPM for all genes,
        % set to 0 to not display anything.


%%%%% OUTPUT
% M2: interaction matrix where the columns corresponding to genes that fall
        % under the expression threshold have been removed.
% intList2: updated list of interacting genes corresponding to the columns of M2
% RNAseq: MATLAB table containing the RNA-seq data for HeLa cells. second
        % column (Gene_Name) holds the gene symbols; fifth column (pTPM) holds the
        % TPM data.

%% load TPM spreadsheet
RNAseq = readtable(RNAseqFileName);

%% plot distribution of TPM is option is selected
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

%% prune columns of M to exclude interactors that fall under the TPM threshold
% better to select non expressed genes because some of the gene names in the
% interactor dataset aren't present in the HeLa gene list so indexing by
% expressed genes would lead to false negatives.
nonExpressedGenes = RNAseq.Gene_name(log(RNAseq.pTPM+1) < expressionThreshold);

expressedGenesIdx = ~ismember(intList,nonExpressedGenes);
intList2 = intList(expressedGenesIdx);
M2 = M(:,expressedGenesIdx);
