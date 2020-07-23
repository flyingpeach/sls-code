function objective = get_L1_norm(transferMtx)
% Induced inf-to-inf objective (i.e. max row sum)

objective = 0;

rowConcat = [];
for k=1:length(transferMtx)
    rowConcat = [rowConcat transferMtx{k}];
end

objective = norm(rowConcat, Inf);
end