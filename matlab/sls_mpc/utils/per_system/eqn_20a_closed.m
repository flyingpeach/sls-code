function [phi_, x_] = eqn_20a_closed(x_loc, psi_, lamb_, y_, z_, ...
                                     cost_, selfIdx, params)
rho = params.rho_;
mu  = params.mu_;

[M1, M2, MSum, MbSum] = eqn_20a_common(x_loc, y_, z_, selfIdx);
a = psi_ - lamb_;
W = pinv(2*((cost_*M2)'*cost_*M2)+rho*(M1'*M1)+mu*MSum)*(rho*M1'*a'+mu*MbSum);

phi_ = (M1*W)';
x_   = M2*W;

end