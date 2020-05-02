function couplingIdx = get_coupling_indices(costMtx)

% TODO: for now, only coupling from cost matrix is considered

couplingIdx = cell(1, length(costMtx));
for i = 1:length(costMtx)
    for j = 1:length(costMtx)
        if costMtx(i,j) ~= 0
            couplingIdx{i} = [couplingIdx{i} j];
        end
    end
end

end