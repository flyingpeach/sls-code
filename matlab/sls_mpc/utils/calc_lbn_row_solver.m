function [v, time] = calc_lbn_row_solver(C, a, B, d)
% Solves general equation 
% min ||v - a||_F^2
%     s.t. B*v + ||C*v||_2 < d
% Tricks used: introduce variable w s.t. w'*w = ||C*v||_2^2, w > 0

dim = length(a); % dimension of v
C1  = [eye(dim) zeros(dim, 1)];
C2  = [zeros(1, dim) 1];

model.Q     = sparse(C1'*C1);
model.obj   = -2*C1'*a;

model.A     = sparse([B 1; zeros(1, dim) -1]);
model.rhs   = [d; 0];
model.sense = '<';

model.quadcon(1).Qc    = sparse(C1'*(C'*C)*C1 - C2'*C2);
model.quadcon(1).q     = sparse(dim+1, 1);
model.quadcon(1).rhs   = 0;

% Note: same effect as setting equality but allows us to use convex solver
model.quadcon(1).sense = '<'; 

LB       = -1e3; % default lb is 0; override
model.lb = LB*ones(size(model.A, 2), 1);

gParams.outputFlag = 0;
result = gurobi(model, gParams);
v      = result.x(1:dim);
time   = result.runtime;

end