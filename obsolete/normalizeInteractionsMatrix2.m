function M = normalizeInteractionsMatrix2(M,normByRows,normByCols)

% sum 
if normByRows
    M = M./repmat(max(M,[],2),1,size(M,2));
end

% z-score transform M by rows and or columns
if normByCols
    M = zscore(M);
end

