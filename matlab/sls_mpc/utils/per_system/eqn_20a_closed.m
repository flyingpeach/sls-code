function [phi_, x_] = eqn_20a_closed(x_loc, psi_, lamb_, y_, z_, ...
                                     cost_, selfIdx, params)
rho = params.rho_;
mu  = params.mu_;
n   = length(x_loc);
nc  = length(z_);    % # of coupled subsystems (including self)

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

W = pinv(2*((cost_*M2)'*cost_*M2)+rho*(M1'*M1)+mu*MSum)*(rho*M1'*a'+mu*MbSum);

phi_ = (M1*W)';
x_   = M2*W;

end