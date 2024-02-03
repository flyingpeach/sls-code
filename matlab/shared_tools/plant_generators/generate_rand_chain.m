function sys = generate_rand_chain(Nx, rho, actDens)
% Populates (A, B2, Nx, Nu) of the specified system with a random chain 
% (tridiagonal A matrix) and a random actuation (B) matrix
% Inputs
%    Nx      : number of states in the system
%    rho     : normalization value; A is generated s.t. max |eig(A)| = rho
%    actDens : actuation density of B, in (0, 1]
%              this is approximate; only exact if things divide exactly
% Outputs
%    sys     : LTISystem containing system matrices

sys    = LTISystem;
sys.Nx = Nx;
sys.Nu = ceil(sys.Nx * actDens);

A                = speye(sys.Nx);
A(1:end-1,2:end) = A(1:end-1,2:end)+ diag(randn(sys.Nx-1, 1));
A(2:end,1:end-1) = A(2:end,1:end-1)+ diag(randn(sys.Nx-1, 1));

B = sparse(sys.Nx, sys.Nu);
for i=1:1:sys.Nu
    B(mod(floor(1/actDens*i-1), sys.Nx)+1,i) = randn();
end

sys.A  = A / max(abs(eigs(A))) * rho;
sys.B2 = B;