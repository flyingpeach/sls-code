function [phi_, time] = mpc_row_solver(x_loc, psi_, lamb_, b1, b2, cost_, rho)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda

n   = length(x_loc);

% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q   = sparse((cost_*x_loc')'*(cost_*x_loc') + rho/2*eye(n));
model.obj = rho*(-psi_ + lamb_);

if isinf(b2)
    model.A     = sparse(x_loc');
    model.rhs   = b1;
    model.sense = '<';
elseif isinf(b1)
    model.A     = sparse(x_loc');
    model.rhs   = b2;
    model.sense = '>';
else
    model.A     = sparse([x_loc'; x_loc']);
    model.rhs   = [b1 b2];
    model.sense = '<>';
end

% default lower bound is 0; override
MPC_LB   = -1e3;
model.lb = MPC_LB*ones(length(model.A), 1);

% solve QP
gParams.outputflag = 0;
result = gurobi(model, gParams);

phi_ = result.x(:);
time = result.runtime;
end