function D = d_step_randomize(M, beta, stabType)
% This is separable for nu, L1, and LInf

n = size(M, 1);

cvx_begin quiet
variable d(n)

if stabType == StabType.L1
    M*d <= beta*d;
    for i=1:n
        d(i) >= 1; % version of d(i) > 1 (numerically easier)
                   % since we can scale D by any constant
    end
elseif stabType == StabType.LInf
    M'*d <= beta*d;
    for i=1:n
        d(i) >= 1;
    end
else % default to nu
    for i=1:n
        for j=1:n
            if M(i,j) ~= 0 && i ~= j
                 log(M(i,j)) + d(i) - d(j) <= log(beta);
            end
        end
    end
end

minimize(0) % arbitrary objective
cvx_end

if stabType == StabType.L1
    D = inv(diag(d));
elseif stabType == StabType.LInf
    D = diag(d);
else % nu
    D = diag(exp(d)); % nu works with log-coordinates, convert back
end
    
end