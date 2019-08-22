function objective = compute_Hinf(R, M, C, D, TFIR)
% return ||[C,D][R;M]||_2->2
% see eqns (3.1), (4.20) of long tutorial

mtx = [];
for t = 1:TFIR
    mtx = blkdiag(mtx, [C,D]*[R{t};M{t}]);
end

objective = sigma_max(mtx);