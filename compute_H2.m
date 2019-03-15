function objective = compute_H2(R,M,C,D,T);
% return ||[C,D][R;M]||_H2^2

objective = 0;
for t = 1:T
    %need to do the vect operation because of quirk in cvx
    vect = vec([C,D]*[R{t};M{t}]);
    objective = objective + vect'*vect;
end