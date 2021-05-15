function [phi_row, time] = mpc_row_solver(x_loc, psi, lamb, ub, lb, cost, rho)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda

n = length(x_loc);

% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q   = sparse((cost*x_loc')'*(cost*x_loc') + rho/2*eye(n));
model.obj = rho*(-psi + lamb);

if isinf(lb)
    model.A     = sparse(x_loc');
    model.rhs   = ub;
    model.sense = '<';
elseif isinf(ub)
    model.A     = sparse(x_loc');
    model.rhs   = lb;
    model.sense = '>';
else
    model.A     = sparse([x_loc'; x_loc']);
    model.rhs   = [ub lb];
    model.sense = '<>';
end

% default lower bound is 0; override
MPC_LB   = -1e3;
model.lb = MPC_LB*ones(length(model.A), 1);

% solve QP
gParams.outputflag = 0;
result = gurobi(model, gParams);

phi_row = result.x(:);
time    = result.runtime;
end