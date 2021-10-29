function phi_row = row_nominal_closed(x_loc, psi, lamb, cost, rho)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda

n   = length(x_loc);
M   = pinv(2*x_loc*cost*cost'*x_loc' + rho*eye(n));

phi_row = rho*(psi - lamb)*M;

end