function objective = get_L1_obj(transferMtx)
objective = 0;
for k=1:length(transferMtx)
    objective = max(objective, norm(transferMtx{k}, Inf));
end
end