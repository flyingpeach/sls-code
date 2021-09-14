function [v, time] = calc_robust_row_solver(C, a, B, d)
% Solves general quadratic equation min ||C*v - a||_F^2 s.t. B*v < d

model.Q     = sparse(C'*C);
model.obj   = -2*C'*a;
model.A     = sparse(B);
model.rhs   = d;
model.sense = '<';

LB       = -1e3; % default lb is 0; override
model.lb = LB*ones(size(model.A, 2), 1);

gParams.outputflag = 0;
result = gurobi(model, gParams);
v      = result.x(:)';
time   = result.runtime;

end