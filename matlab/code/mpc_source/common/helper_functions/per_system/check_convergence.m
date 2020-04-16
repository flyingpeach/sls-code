function converged = check_convergence(phi_loc, psi_loc, psi_prev_loc, eps_p, eps_d)
    
    primDiff  = norm(phi_loc-psi_loc,'fro');
    dualDiff  = norm(psi_loc-psi_prev_loc,'fro');
    converged = primDiff <= eps_p && dualDiff <= eps_d;

end
