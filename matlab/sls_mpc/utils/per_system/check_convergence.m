function [converged, scale] = check_convergence(prim1, prim2, dual1, dual2, params)
% Checks whether the following is satisfied, where ||  ||_F is frob norm
% || prim1 - prim2 ||_F <= params.eps_p_
% || dual1 - dual2 ||_F <= params.eps_d_
% Returns scaling factor for adaptive ADMM rho
% If adaptive ADMM is not specified, scale = 1 (i.e. rho does not change)

primRes   = norm(prim1 - prim2, 'fro');
dualRes   = norm(dual1 - dual2, 'fro');
converged = primRes <= params.eps_p_ && dualRes <= params.eps_d_;

scale = 1;
if params.has_adaptive_admm()
     if primRes > params.muAdapt_ * dualRes
         scale = params.tau_i_;
     elseif dualRes > params.muAdapt_ * primRes
         scale = 1 / params.tau_d_;
     end
end

end