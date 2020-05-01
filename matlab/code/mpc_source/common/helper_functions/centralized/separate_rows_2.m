function [Psi_rows, Lambda_rows] = separate_rows_2(sys, tFIR, ...
                                                        r, s_r, r_loc, m_loc, ...
                                                        Psi, Lambda)
Nx = sys.Nx; Nu = sys.Nu;

Psi_rows    = cell(1, Nx);
Lambda_rows = cell(1, Nx);
                
% Separate the given matrices
k = 0;
for i_ = 1:Nx
    if mod(i_, Nx/Nu) == 0
        k = k+1;
    end
    
    for i = r{i_}
        j = find(r{i_}==i);
        if j<=tFIR
            Psi_rows{i} = Psi(i,s_r{i_}(j,1:max(length(find(r_loc(i_,:))))));
            Lambda_rows{i} = Lambda(i,s_r{i_}(j,1:max(length(find(r_loc(i_,:))))));
        else
            Psi_rows{i} = Psi(i,s_r{i_}(j,1:max(length(find(m_loc(k,:))))));
            Lambda_rows{i} = Lambda(i,s_r{i_}(j,1:max(length(find(m_loc(k,:))))));
        end
    end
end

end