function [Phi_cols, Lambda_cols] = separate_cols(sys, c, s_c, Phi, Lambda)

Nx = sys.Nx;

Phi_cols    = cell(Nx, 1);
Lambda_cols = cell(Nx, 1);               

for i = 1:Nx
    Phi_cols{i}    = Phi(s_c{i}, c{i});
    Lambda_cols{i} = Lambda(s_c{i}, c{i});
end
        
