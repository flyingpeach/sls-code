function [phi_, x_] = mpc_coupled_row_closed(x_loc, psi, lamb, y, z, ...
                                             cost, sIdx, params)
rho = params.rho_;
mu  = params.mu_;

[M1, M2, MSum, MbSum] = coupled_row_setup(x_loc, y, z, sIdx);
a = psi - lamb;
W = pinv(2*((cost*M2)'*cost*M2)+rho*(M1'*M1)+mu*MSum)*(rho*M1'*a'+mu*MbSum);

phi_ = (M1*W)';
x_   = M2*W;

end