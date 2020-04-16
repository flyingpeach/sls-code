function Phi_loc = eqn_16a_closed(s_ri, x_ri, psi_rowi, lamb_rowi, rho)

M       = inv(2 * x_ri * x_ri' + rho * eye(size(s_ri, 2)));
Phi_loc = rho * (psi_rowi - lamb_rowi) * M;

end