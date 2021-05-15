function converged = check_convergence_cons(z_cp, x, z, z_prev, params)

primRes   = norm(x - z_cp, 'fro');
dualRes   = norm(z - z_prev, 'fro');
converged = primRes <= params.eps_x_ && dualRes <= params.eps_z_;

end
