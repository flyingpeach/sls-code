function objective = get_eps1_obj(transferMtx)
% Element-wise l1 norm
objective = 0;
for k=1:length(transferMtx)
    objective = objective + sum(vec(abs(transferMtx{k})));
end

end