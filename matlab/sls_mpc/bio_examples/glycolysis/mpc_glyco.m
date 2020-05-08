function [x_nxt, u, status] = mpc_glyco(sys, tFIR, x_lb, u_lb, u_ub, xt, relax)
% Note that these are "x, u" in the local sense of mpc
% In the bilophila example, they are y and u_tilde (shifted coordinates)

Nx = sys.Nx; Nu = sys.Nu; A = sys.A; B = sys.B2;

RLoc = ones(Nx, Nx); % no sparsity constraint
MLoc = ones(Nu, Nx);

count = 0;
for k = 1:tFIR
    RSupp{k} = RLoc;
    MSupp{k} = MLoc;
    count = count + sum(sum(RSupp{k})) + sum(sum(MSupp{k}));
end

cvx_begin

variable X(count)
expression Rs(Nx, Nx, tFIR)
expression Ms(Nu, Nx, tFIR)
   
R = cell(tFIR, 1);
M = cell(tFIR, 1);

for t = 1:tFIR
    R{t} = Rs(:,:,t); M{t} = Ms(:,:,t); 
end

[R, M] = add_sparse_constraints(R, M, RSupp, MSupp, X, tFIR);

objective = 0;

input_pen = 1e2;
for k = 1:tFIR
      objective = objective + [0 0 1]*R{k}*xt;
      objective = objective - input_pen*[1 1 1]*M{k}*xt; 
      % since u is always positive
end

maximize(objective)
subject to

% Achievability constraints
R{1} == eye(Nx); 
for k=1:tFIR-1
    R{k+1} == A*R{k} + B*M{k};
end

for k=1:tFIR
     M{k}*xt >= u_lb;
     M{k}*xt <= u_ub;
end

kstart=1;
if relax
    % bypass R{1}*xt if it causes infeasibility
    kstart=2;
end
for k=kstart:tFIR
     R{k}*xt >= x_lb;
end

cvx_end

u     = M{1}*xt;
x_nxt = R{2}*xt;

status = cvx_status;

end
