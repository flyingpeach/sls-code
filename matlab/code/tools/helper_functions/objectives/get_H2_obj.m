function objective = get_H2_obj(sys, R, M)
objective = 0;
for t = 1:length(R)
    % need to do the vect operation because of quirk in cvx
    vect = vec([sys.C1, sys.D12]*[R{t};M{t}]);
    objective = objective + vect'*vect;
end
end
