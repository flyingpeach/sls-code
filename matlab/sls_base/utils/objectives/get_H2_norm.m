function objective = get_H2_norm(transferMtx)
% Frobenius
objective = 0;
for k=1:length(transferMtx)
    % need to do the vect operation because of quirk in cvx
    vect = vec(transferMtx{k});
    objective = objective + vect'*vect;
end
end
