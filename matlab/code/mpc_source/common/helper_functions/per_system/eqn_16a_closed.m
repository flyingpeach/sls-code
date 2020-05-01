function Phi_loc = eqn_16a_closed(x_ri, psi_rowi, lamb_rowi, params)

n   = length(x_ri);
rho = params.rho_;

M       = inv(2 * x_ri * x_ri' + rho * eye(n));
Phi_loc = rho * (psi_rowi - lamb_rowi) * M;

end