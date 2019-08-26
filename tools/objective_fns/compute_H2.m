function objective = compute_H2(sys, params, R, M)
% return ||[C1,D12][R;M]||_H2^2 as per (4.20)

objective = 0;
for t = 1:params.tFIR_
    %need to do the vect operation because of quirk in cvx
    vect = vec([sys.C1, sys.D12]*[R{t};M{t}]);
    objective = objective + vect'*vect;
end