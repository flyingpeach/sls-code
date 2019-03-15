function [A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,act_dens);
% Nx = size of chain
% alpha is how much state is spread between neighbors
% rho is stabiilty of A.  Choose rho = 1 for dbl-stoch A
% act_dens \in (0,1] is atuation density.  This is approximate (only exact
% if things divide into each other evenly)
% x_1(t+1) = rho*[(1-alpha)*x_1(t) + alpha x_2(t)] + B(1,1)u_1(t)
% x_i(t+1) = rho*[alpha*x_{i-1}(t) + (1-2*alpha)x_i(t) + alpha*x_{i+1}(t)]
% + B(i,i)u_i(t)
% x_N(t+1) = rho*[alpha*x_{N-1}(t) + (1-alpha)x_N(t)] + B(N,N)u_N(t)


if (nargin < 1)
    alpha = .2;
    rho = 1;
    Nx = 9;
    act_dens = 1/3;
end

Nu = ceil(Nx*act_dens);

A = (1-2*alpha)*speye(Nx);
A(1,1) = A(1,1) + alpha;
A(Nx,Nx) = A(Nx,Nx) + alpha;
A(1:end-1,2:end) = A(1:end-1,2:end)+ alpha*speye(Nx-1);
A(2:end,1:end-1) = A(2:end,1:end-1)+ alpha*speye(Nx-1);

A = A*rho;

B = sparse(Nx, Nu);
for i=1:1:Nu
    B(mod(floor(1/act_dens*i-1),Nx)+1,i) = 1;
end
