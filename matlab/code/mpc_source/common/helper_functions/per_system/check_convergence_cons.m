function converged = check_convergence_cons(z_cp, x_, z_, z_prev_, params)

eps_x = params.eps_x_;
eps_z = params.eps_z_;

primRes   = norm(x_-z_cp, 'fro');
dualRes   = norm(z_-z_prev_, 'fro');
converged = primRes <= eps_x && dualRes <= eps_z;

end
