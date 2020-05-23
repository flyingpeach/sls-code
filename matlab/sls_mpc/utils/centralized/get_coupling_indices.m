function cpIdx = get_coupling_indices(C, K)
% C is cost matrix
% K is constraint matrix
% Both were augmented from state/input matrices to be the same size as Phi

nVals = length(C);
cpIdx = cell(1, nVals);

for i = 1:nVals
    for j = 1:nVals
        if (C(i,j) ~= 0 || K(i,j) ~=0)
            cpIdx{i} = unique([cpIdx{i} j]);
            cpIdx{j} = unique([cpIdx{j} i]);
        end
    end
end

end