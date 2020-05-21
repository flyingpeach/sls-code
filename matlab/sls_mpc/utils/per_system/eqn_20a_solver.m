function [phi_, x_] = eqn_20a_solver(x_loc, psi_, lamb_, y_, z_, ...
                                     cost_, k_, selfIdx, params)        
rho = params.rho_;
mu  = params.mu_;
n   = length(x_loc);
nc  = length(z_);

MSum  = 0; 
MbSum = 0;
M     = cell(nc, 1);
for j = 1:nc
    if j < selfIdx
            M{j}        = zeros(1,nc+n-1); 
            M{j}(n+j)   = 1;
    elseif j == selfIdx
            M{j}        = [x_loc' zeros(1,nc-1)];
    elseif j > selfIdx
            M{j}        = zeros(1,nc+n-1); 
            M{j}(n+j-1) = 1;
    end

    MSum  = MSum + (M{j}'*M{j});

    b     = z_{j} - y_{j};    
    MbSum = MbSum + M{j}'*b;
end

a     = psi_ - lamb_;
myEye = eye(nc-1);
M1    = [eye(n) zeros(n,nc-1)];
M2    = [zeros(selfIdx-1,n)  myEye(1:selfIdx-1,:); 
         x_loc'            zeros(1,nc-1); 
         zeros(nc-selfIdx,n) myEye(selfIdx:end,:)];

% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q   = sparse((cost_*M2)'*(cost_*M2) + rho/2*(M1'*M1) + mu/2*(MSum));
model.obj = -rho*a*M1 -mu*MbSum';

ub  = params.stateUB_;
lb  = params.stateLB_;

model.A     = sparse(k_*M2);
if isinf(lb)
    model.rhs   = [ub ub];
    model.sense = '<';
elseif isinf(ub)
    model.rhs   = [lb lb];
    model.sense = '>';
else
    model.rhs   = [ub lb]; 
    model.sense = '<>';
end

% default lower bound is 0; override
MPC_LB   = -100;
model.lb = MPC_LB*ones(nc+n-1,1);

% solve QP
gParams.outputflag = 0;
result = gurobi(model, gParams);
W      = result.x(:); 

phi_  = (M1*W)';
x_    = M2*W;
end