function [Phi_cols, Lambda_cols] = separate_cols(c, s_c, Nx)

Phi_cols    = cell(1, Nx);
Lambda_cols = cell(1, Nx);               

for i = 1:Nx
    Phi_cols{i}    = Phi(s_c{i}, c{i});
    Lambda_cols{i} = Lambda(s_c{i}, c{i});
end
        
