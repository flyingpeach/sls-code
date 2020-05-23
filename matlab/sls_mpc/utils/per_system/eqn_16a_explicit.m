function phi_ = eqn_16a_explicit(x_loc, psi_, lamb_, b1, b2, cost_, rho)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda

n   = length(x_loc);
M   = inv(2*x_loc*cost_*cost_'*x_loc' + rho*eye(n));
a   = psi_ - lamb_;
    
crit  = rho*a*M*x_loc;

lamb = 0;
if crit - b1 > 0
    lamb = (crit - b1) / (x_loc'*M*x_loc);
elseif crit - b2 < 0
    lamb = (crit - b2) / (x_loc'*M*x_loc);
end

phi_ = (rho*a - lamb*x_loc')*M;

end