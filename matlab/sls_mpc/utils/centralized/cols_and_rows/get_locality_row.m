function s_r = get_locality_row(SuppMtx)
% s_r{i} gives the set of nonzero columns associated with ith row

nRows = size(SuppMtx, 1);
s_r   = cell(nRows, 1);

for row=1:nRows
    s_r{row} = find(SuppMtx(row, :));
end

end