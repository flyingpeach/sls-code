function objective = get_1to1_norm(transferMtx)
% Induced 1-to-1 objective (i.e. max column sum)
objective = 0;
for k=1:length(transferMtx)
    objective = max(objective, norm(transferMtx{k}, 1));
end
end