function [Phi_cols, Lambda_cols] = separate_cols(sys, c, s_c, Phi, Lambda)

Nx = sys.Nx; Nu = sys.Nu;

Phi_cols    = cell(1, Nx);
Lambda_cols = cell(1, Nx);               

for i = 1:Nx
    Phi_cols{i}    = Phi(s_c{i}, c{i});
    Lambda_cols{i} = Lambda(s_c{i}, c{i});
end
        
