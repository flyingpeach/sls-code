function converged = check_convergence_cons(z_av, x_loc, z_loc, z_prev_loc, eps_x, eps_z)
    
    primRes   = norm(x_loc-z_av, 'fro');
    dualRes   = norm(z_loc-z_prev_loc, 'fro');
    converged = primRes <= eps_x && dualRes <= eps_z;

end
