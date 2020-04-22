function [Phi_loc, X_loc] = eqn_20a_closed(x_ri, psi_rowi, lamb_rowi, ...
                                           Z_locs, y_rowi, ...
                                           indicesi, i_new, ci, n, rho, mu)

m_j   = max(length(indicesi));
a     = psi_rowi - lamb_rowi;

Eye = eye(m_j-1);
M1  = [eye(n) zeros(n,m_j-1)];
M2  = [zeros(i_new-1,n) Eye(1:i_new-1,:); x_ri' zeros(1,m_j-1); zeros(m_j-i_new,n) Eye(i_new:end,:)];

Mj_sum  = 0; 
Mjb_sum = 0;
M       = cell(1, m_j);
for j = 1:m_j
    if j < i_new
            M{j}      = zeros(1,m_j+n-1); 
            M{j}(n+j) = 1;
    elseif j == i_new
            M{j} = [x_ri' zeros(1,m_j-1)];
    elseif j>i_new
            M{j}        = zeros(1,m_j+n-1); 
            M{j}(n+j-1) = 1;
    end

    Mj_sum = Mj_sum + (M{j}'*M{j});
    
    k       = indicesi(j);
    b       = Z_locs{k} - y_rowi{k};    
    Mjb_sum = Mjb_sum + M{j}'*b;
end

W = pinv(2*((ci*M2)'*ci*M2)+rho*(M1'*M1)+mu*Mj_sum)*(rho*M1'*a'+mu*Mjb_sum);

Phi_loc         = (M1*W)';
X_loc           = zeros(Nx*tFIR+Nu*(tFIR-1),1); 
X_loc(indicesi) = M2*W;

end