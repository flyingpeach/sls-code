function [phi_, x_] = eqn_20a_solver(x_loc, psi_, lamb_, y_, z_, ...
                                     cost_, k_, selfIdx, lb, ub, params)        
rho = params.rho_;
mu  = params.mu_;

[M1, M2, MSum, MbSum] = eqn_20a_common(x_loc, y_, z_, selfIdx);
a = psi_ - lamb_;

% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q   = sparse((cost_*M2)'*(cost_*M2) + rho/2*(M1'*M1) + mu/2*(MSum));
model.obj = -rho*a*M1 -mu*MbSum';

if isinf(lb)
    model.A     = sparse(k_*M2);
    model.rhs   = ub;
    model.sense = '<';
elseif isinf(ub)
    model.A     = sparse(k_*M2);
    model.rhs   = lb;
    model.sense = '>';
else
    model.A     = sparse([k_; k_]*M2);
    model.rhs   = [ub lb]; 
    model.sense = '<>';
end

% default lower bound is 0; override
MPC_LB   = -100;
model.lb = MPC_LB*ones(length(model.A), 1);

% solve QP
gParams.outputflag = 0;
result = gurobi(model, gParams);
W      = result.x(:); 

phi_  = (M1*W)';
x_    = M2*W;
end