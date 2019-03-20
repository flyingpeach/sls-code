function [A,B] = generate_rand_chain(Nx,rho,act_dens);
% Nx = size of chain
% generates a random chain (tridiagonal A matrix) with random B matrix of
% specified actuation density. Normlaized so max |eig(A)| = rho;

Nu = ceil(Nx*act_dens);

A = speye(Nx);
A(1:end-1,2:end) = A(1:end-1,2:end)+ diag(randn(Nx-1));
A(2:end,1:end-1) = A(2:end,1:end-1)+ diag(randn(Nx-1));

A = A/max(abs(eigs(A)))*rho;

B = sparse(Nx, Nu);
for i=1:1:Nu
    B(mod(floor(1/act_dens*i-1),Nx)+1,i) = randn();
end
