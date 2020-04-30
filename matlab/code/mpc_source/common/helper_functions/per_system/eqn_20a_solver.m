function [Phi_loc, X_loc] = eqn_20a_solver(x_ri, psi_rowi, lamb_rowi, ...
                                                 Z_locs, y_rowi, ki, ...
                                                 indicesi, i_new, ci, n, rho, mu, up)

m_j = max(length(indicesi));
a   = psi_rowi - lamb_rowi;

myEye = eye(m_j-1);
M1  = [eye(n) zeros(n,m_j-1)];
M2  = [zeros(i_new-1,n) myEye(1:i_new-1,:); x_ri' zeros(1,m_j-1); zeros(m_j-i_new,n) myEye(i_new:end,:)];

Mj_sum  = 0; 
Mjb_sum = 0;
M       = cell(1, m_j);
for j = 1:m_j
    if j < i_new
            M{j}      = zeros(1,m_j+n-1); 
            M{j}(n+j) = 1;
    elseif j == i_new
            M{j} = [x_ri' zeros(1,m_j-1)];
    elseif j > i_new
            M{j}        = zeros(1,m_j+n-1); 
            M{j}(n+j-1) = 1;
    end

    Mj_sum = Mj_sum + (M{j}'*M{j});
    
    k       = indicesi(j);
    b       = Z_locs{k} - y_rowi{k};    
    Mjb_sum = Mjb_sum + M{j}'*b;
end

% set up QP
% minimize Phi*Q*Phi' + obj*Phi'
% note: constant terms omitted since we are interested in argmin only
model.Q   = sparse((ci*M2)'*(ci*M2) + rho/2*(M1'*M1) + mu/2*(Mj_sum));
model.obj = -rho*a*M1 -mu*Mjb_sum';

model.A     = sparse(ki*M2);
model.rhs   = [up up]; 
model.sense = '<';

% default lower bound is 0; override
MPC_LB   = -100;
model.lb = MPC_LB*ones(m_j+n-1,1);

% solve QP
gParams.outputflag = 0;
result = gurobi(model, gParams);
W      = result.x(:); 

Phi_loc         = (M1*W)';
X_loc           = zeros(length(Z_locs), 1); 
X_loc(indicesi) = M2*W;
end