function objective = get_stab_norm(transferMtx)
% Stability objective used to produce plots in previous publications
bigMtx = [];

for k=1:length(transferMtx)
    bigMtx = [bigMtx, transferMtx{k}];
end

objective = norm(bigMtx, Inf);

end