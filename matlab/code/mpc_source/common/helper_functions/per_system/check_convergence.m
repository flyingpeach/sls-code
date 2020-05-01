function converged = check_convergence(phi_, psi_, psi_prev_, params)
% phi_      : rows of Phi
% psi_      : rows of Psi
% psi_prev_ : rows of Psi from previous ADMM iteration
% params    : MPC parameters

eps_d = params.eps_d_;
eps_p = params.eps_p_;

primRes   = norm(phi_-psi_,'fro');
dualRes   = norm(psi_-psi_prev_,'fro');
converged = primRes <= eps_p && dualRes <= eps_d;

end
