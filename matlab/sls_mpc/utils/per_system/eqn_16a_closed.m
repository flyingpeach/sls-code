function phi_ = eqn_16a_closed(x_loc, psi_, lamb_, rho)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda

n   = length(x_loc);
M   = inv(2*x_loc*x_loc' + rho*eye(n));

phi_ = rho*(psi_ - lamb_)*M;

end