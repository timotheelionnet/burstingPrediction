function M = normalizeInteractionsMatrix(M,normByRows,normByCols)

% z-score transform M by rows and or columns
if normByCols
    M = zscore(M);
end

if normByRows
    M = zscore(M')';
end



