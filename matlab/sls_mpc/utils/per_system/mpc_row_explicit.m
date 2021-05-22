function phi_row = mpc_row_explicit(x_loc, psi, lamb, ub, lb, cost, rho)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda

n   = length(x_loc);
M   = inv(2*x_loc*cost*cost'*x_loc' + rho*eye(n));
a   = psi - lamb;
    
crit  = rho*a*M*x_loc;

lamb = 0;
if crit - ub > 0
    lamb = (crit - ub) / (x_loc'*M*x_loc);
elseif crit - lb < 0
    lamb = (crit - lb) / (x_loc'*M*x_loc);
end

phi_row = (rho*a - lamb*x_loc')*M;

end