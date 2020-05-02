function [phi_, x_] = eqn_20a_solver(x_loc, psi_, lamb_, y_, Z_rows, ...
                                     k_, c_, cp_, i_, params)

n   = length(x_loc);
rho = params.rho_;
mu  = params.mu_;
up  = params.constrUpperbnd_;

m_j = max(length(cp_));
a   = psi_ - lamb_;

myEye = eye(m_j-1);
M1  = [eye(n) zeros(n,m_j-1)];
M2  = [zeros(i_-1,n) myEye(1:i_-1,:); x_loc' zeros(1,m_j-1); zeros(m_j-i_,n) myEye(i_:end,:)];

Mj_sum  = 0; 
Mjb_sum = 0;
M       = cell(1, m_j);
for j = 1:m_j
    if j < i_
            M{j}      = zeros(1,m_j+n-1); 
            M{j}(n+j) = 1;
    elseif j == i_
            M{j} = [x_loc' zeros(1,m_j-1)];
    elseif j > i_
            M{j}        = zeros(1,m_j+n-1); 
            M{j}(n+j-1) = 1;
    end

    Mj_sum = Mj_sum + (M{j}'*M{j});
    
    k       = cp_(j);
    b       = Z_rows{k} - y_{k};    
    Mjb_sum = Mjb_sum + M{j}'*b;
end

% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q   = sparse((c_*M2)'*(c_*M2) + rho/2*(M1'*M1) + mu/2*(Mj_sum));
model.obj = -rho*a*M1 -mu*Mjb_sum';

model.A     = sparse(k_*M2);
model.rhs   = [up up]; 
model.sense = '<';

% default lower bound is 0; override
MPC_LB   = -100;
model.lb = MPC_LB*ones(m_j+n-1,1);

% solve QP
gParams.outputflag = 0;
result = gurobi(model, gParams);
W      = result.x(:); 

phi_         = (M1*W)';
x_           = zeros(length(Z_rows), 1); 
x_(cp_) = M2*W;
end