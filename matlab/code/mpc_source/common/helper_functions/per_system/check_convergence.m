function converged = check_convergence(phi_loc, psi_loc, psi_prev_loc, eps_p, eps_d)
    
    primRes   = norm(phi_loc-psi_loc,'fro');
    dualRes   = norm(psi_loc-psi_prev_loc,'fro');
    converged = primRes <= eps_p && dualRes <= eps_d;

end
