function objective = get_HInf_obj(transferMtx)
bigMtx = [];
for k=1:length(transferMtx)
    bigMtx = blkdiag(bigMtx, transferMtx{k});
end

objective = sigma_max(full(bigMtx));
end