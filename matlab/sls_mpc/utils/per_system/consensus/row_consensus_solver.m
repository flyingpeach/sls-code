function [phi_row, x_, time] = row_consensus_solver(x_loc, psi, lamb, y, z, ...
                                   cost, constr, sIdx, ub, lb, params)        
rho = params.rho_;
mu  = params.mu_;

[M1, M2, MSum, MbSum] = row_consensus_setup(x_loc, y, z, sIdx);
a = psi - lamb;

% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q   = sparse((cost*M2)'*(cost*M2) + rho/2*(M1'*M1) + mu/2*(MSum));
model.obj = -rho*a*M1 -mu*MbSum';

if isinf(ub)
    model.A     = sparse(constr*M2);
    model.rhs   = lb;
    model.sense = '<';
elseif isinf(lb)
    model.A     = sparse(constr*M2);
    model.rhs   = ub;
    model.sense = '>';
else
    model.A     = sparse([constr; constr]*M2);
    model.rhs   = [lb ub]; 
    model.sense = '<>';
end

% default lower bound is 0; override
MPC_LB   = -1e3;
model.lb = MPC_LB*ones(length(model.A), 1);

% solve QP
gParams.outputflag = 0;
result = gurobi(model, gParams);
W      = result.x(:); 

phi_row  = (M1*W)';
x_    = M2*W;
time  = result.runtime;
end