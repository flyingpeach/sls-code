function objective = get_HInf_obj(sys, R, M)
mtx = [];
for t = 1:length(R)
    mtx = blkdiag(mtx, [sys.C1, sys.D12]*[R{t};M{t}]);
end

objective = sigma_max(full(mtx));
end