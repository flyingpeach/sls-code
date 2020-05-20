function [phi_, x_] = eqn_20a_closed(x_loc, psi_, lamb_, y_, Z_rows, ...
                                     cost_, coupl_, self_, params)
             
rho = params.rho_;
mu  = params.mu_;                                   

n   = length(x_loc);
nc  = length(coupl_); % number coupling

MSum  = 0; 
MbSum = 0;
M     = cell(nc, 1);
for j = 1:nc
    if j < self_
            M{j}        = zeros(1,nc+n-1); 
            M{j}(n+j)   = 1;
    elseif j == self_
            M{j}        = [x_loc' zeros(1,nc-1)];
    elseif j > self_
            M{j}        = zeros(1,nc+n-1); 
            M{j}(n+j-1) = 1;
    end

    MSum  = MSum + (M{j}'*M{j});

    k     = coupl_(j);
    b     = Z_rows{k} - y_{k};    
    MbSum = MbSum + M{j}'*b;
end

a     = psi_ - lamb_;
myEye = eye(nc-1);
M1    = [eye(n) zeros(n,nc-1)];
M2    = [zeros(self_-1,n)  myEye(1:self_-1,:); 
         x_loc'            zeros(1,nc-1); 
         zeros(nc-self_,n) myEye(self_:end,:)];

W = pinv(2*((cost_*M2)'*cost_*M2)+rho*(M1'*M1)+mu*MSum)*(rho*M1'*a'+mu*MbSum);

phi_       = (M1*W)';
x_         = zeros(length(Z_rows), 1); 
x_(coupl_) = M2*W;

end