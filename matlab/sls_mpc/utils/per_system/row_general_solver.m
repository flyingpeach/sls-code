function [v, time] = row_general_solver(C, a, F, B, d)
% Solves general quadratic equation min ||C*v - a||^2 + F*v s.t. B*v < d
% Problem formulation: the last variable corresponds to consensus 
% variable, so return that as v2 (and the rest as v1)

model.Q     = sparse(C'*C);
model.obj   = F' - 2*C'*a;
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