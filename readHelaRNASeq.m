function [M,intList,RNAseq] = readHelaRNASeq(RNAseqFileName,intList,expressionThreshold, M)

RNAseq = readtable(RNAseqFileName);
%%
%isInIntList = ismember(RNAseq.Gene_name,intList);

%rnaLevels = RNAseq.pTPM(isInIntList);

[n1,x1] = hist(log(RNAseq.pTPM+1),0:.1:10);
figure('Name',['Expression levels in Hela, threshold set to ',num2str(expressionThreshold)]);
hold;
plot(x1,n1/sum(n1),'DisplayName','All genes');
[n2,x2] = hist(log(RNAseq.pTPM(ismember(RNAseq.Gene_name,intList))+1),0:.1:10);
plot(x2,n2/sum(n2),'DisplayName','Interactors');
xlabel('log(pTPM+1)');
ylabel('Counts');

%%
% better to select non expressed genes because some of the gene names in the
% interactor dataset aren't present in the HeLa gene list so indexing by
% expressed genes would lead to false negatives.
nonExpressedGenes = RNAseq.Gene_name(log(RNAseq.pTPM+1) < expressionThreshold);

expressedGenesIdx = ~ismember(intList,nonExpressedGenes);
intList = intList(expressedGenesIdx);
M = M(:,expressedGenesIdx);
