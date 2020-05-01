function phi_ = eqn_16a_solver(x_loc, psi_, lamb_, params)
% x_loc  : locally observed state
% psi_   : row of Psi
% lamb_  : row of Lambda
% params : MPC parameters

n   = length(x_loc);
rho = params.rho_;
    
% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q     = sparse(x_loc*x_loc' + rho/2*eye(n));
model.obj   = rho*(-psi_ + lamb_);

% s.t. x_ri'*Phi < upperbnd
%      x_ri'*Phi > lowerbnd
model.A     = sparse([x_loc'; x_loc']);
model.rhs   = [params.stateUpperbnd_, params.stateLowerbnd_];
model.sense = '<>';
    
% default lower bound is 0; override
MPC_LB   = -100;
model.lb = MPC_LB*ones(n, 1);
    
% solve QP    
gParams.outputflag = 0;
result = gurobi(model, gParams);
phi_   = result.x(:);
end