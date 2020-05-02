function objective = get_L1_norm(transferMtx)
% Induced inf-to-inf objective (i.e. max row sum)
objective = 0;
for k=1:length(transferMtx)
    objective = max(objective, norm(transferMtx{k}, Inf));
end
end