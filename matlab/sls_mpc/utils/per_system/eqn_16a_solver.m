function phi_ = eqn_16a_solver(x_loc, psi_, lamb_, b1, b2, cost_, rho)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda

cvx_begin quiet
variable phi(size(psi_))

objVec = vec(cost_ * phi * x_loc);
admmVec = vec(phi - psi_ + lamb_);

objective = objVec'*objVec + 0.5*rho*admmVec'*admmVec;

minimize objective
subject to

if ~isinf(b1)
    phi*x_loc <= b1;
end
if ~isinf(b2)
    phi*x_loc >= b2;
end

cvx_end

phi_ = phi; % cvx doesn't allow variable with underscore ending

end