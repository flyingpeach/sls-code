function [D, beta] = d_step_minimize(M, stabType)
% Note: this is only separable for nu
% We include non-separable formulations for L1, LInf to compare

n = size(M, 1);

if ismember(stabType, [StabType.L1, StabType.LInf])
    if stabType == StabType.L1
        [eigvecs, eigvals] = eig(M);
    else % LInf
        [eigvecs, eigvals] = eig(M');
    end
    eigvals = eigvals(eigvals ~= 0); % convert to vector   
    specRad = max(abs(eigvals));
    % locations of max magnitude (right, left) eigenvalue(s)
    maxIdxs = find(abs(eigvals) == specRad); 
    
    eigvec = eigvecs(:,maxIdxs(1));
    
    if stabType == StabType.L1
        D = inv(diag(eigvec));
    else % LInf
        D = diag(eigvec);
    end
            
else % default to nu
    cvx_begin quiet
    variable eta 
    variable l(n)

    for i=1:n
        for j=1:n
            if M(i,j) ~= 0
                if i ~= j
                    log(M(i,j)) + l(i) - l(j) <= eta;
                else % i == j
                    log(M(i,j)) <= eta;
                end
            end
        end
    end

    minimize(eta)
    cvx_end
    D = diag(exp(l)); % nu works with log-coordinates, convert back
end

beta = get_bound(D, M, stabType);

end