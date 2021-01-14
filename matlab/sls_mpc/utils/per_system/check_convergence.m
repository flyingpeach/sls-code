function converged = check_convergence(prim1, prim2, dual1, dual2, params)
% Checks whether the following is satisfied, where ||  ||_F is frob norm
% || prim1 - prim2 ||_F <= params.eps_p_
% || dual1 - dual2 ||_F <= params.eps_d_

primRes = norm(prim1 - prim2, 'fro');
dualRes = norm(dual1 - dual2, 'fro');
converged = primRes <= params.eps_p_ && dualRes <= params.eps_d_;

end