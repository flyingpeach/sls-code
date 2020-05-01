function [phi_, x_] = eqn_20a_closed(x_loc, psi_, lamb_, y_, Z_rows, ...
                                           c_, cp_, i_, params)

n   = length(x_loc);
rho = params.rho_;
mu  = params.mu_;                                       
                                       
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

W = pinv(2*((c_*M2)'*c_*M2)+rho*(M1'*M1)+mu*Mj_sum)*(rho*M1'*a'+mu*Mjb_sum);

phi_         = (M1*W)';

x_           = zeros(length(Z_rows), 1); 
x_(cp_) = M2*W;

end