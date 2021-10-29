function s_c = get_locality_col(SuppMtx)
% s_c{i} gives the set of nonzero rows associated with the ith column

nCols = size(SuppMtx, 2);
s_c   = cell(nCols, 1);

for col=1:nCols
    s_c{col} = find(SuppMtx(:, col));
end

end