function objective = compute_Hinf(sys, params, R, M)
% return max singular value of [C1,D12][R;M]

mtx = [];
for t = 1:params.tFIR_
    mtx = blkdiag(mtx, [sys.C1, sys.D12]*[R{t};M{t}]);
end

objective = sigma_max(mtx);