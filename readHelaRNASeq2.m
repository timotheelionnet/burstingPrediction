function [tf2,RNAseq] = readHelaRNASeq2(RNAseqFileName,tf,expressionThreshold)

RNAseq = readtable(RNAseqFileName);
%%

%%
% better to select non expressed genes because some of the gene names in the
% interactor dataset aren't present in the HeLa gene list so indexing by
% expressed genes would lead to false negatives.
nonExpressedGenes = RNAseq.Gene_name(log(RNAseq.pTPM+1) < expressionThreshold);
tf2 = tf;
tf2.expressedInHela = setdiff(tf.int,nonExpressedGenes);
tf2.trainInt = intersect(tf.trainInt,tf2.expressedInHela);
