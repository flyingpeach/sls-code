function objective = get_L1_obj(sys, R, M)
mtx = [];
for t = 1:length(R)
    mtx = blkdiag(mtx, [sys.C1, sys.D12]*[R{t};M{t}]);
end

objective = norm(mtx, Inf); % note: L1 is induced inf-inf norm
end