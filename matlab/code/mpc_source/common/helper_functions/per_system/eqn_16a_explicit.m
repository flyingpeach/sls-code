function phi_ = eqn_16a_explicit(x_loc, psi_, lamb_, params)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda
% params : MPC parameters

n   = length(x_loc);
rho = params.rho_;
M   = inv(2*x_loc*x_loc' + rho*eye(n));
a   = psi_ - lamb_;
    
b1 = params.stateUpperbnd_; 
b2 = params.stateLowerbnd_;

crit  = rho*a*M*x_loc;

lamb = 0;
if crit - b1 > 0
    lamb = (crit - b1) / (x_loc'*M*x_loc);
elseif crit - b2 < 0
    lamb = (crit - b2) / (x_loc'*M*x_loc);
end

phi_ = (rho*a - lamb*x_loc')*M;

end