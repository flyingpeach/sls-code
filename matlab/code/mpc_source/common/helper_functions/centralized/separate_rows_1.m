function [Psi_rows, Lambda_rows] = separate_rows_1(sys, tFIR, ...
                                                        r, s_r, r_loc, m_loc, ...
                                                        Psi, Lambda)
Nx = sys.Nx; Nu = sys.Nu;

Psi_rows    = cell(1, Nx);
Lambda_rows = cell(1, Nx);
        
k = 0;
for i = 1:Nx
    if mod(i, Nx/Nu) == 0
        k = k+1;
        for j = 1:tFIR+(tFIR-1)
            if j<=tFIR
                Psi_rows{i}    = Psi(r{i},s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))));
                Lambda_rows{i} = Lambda(r{i},s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))));
            else
                Psi_rows{i}    = Psi(r{i},s_r{i}(j,1:max(length(find(m_loc{j-tFIR}(k,:))))));
                Lambda_rows{i} = Lambda(r{i},s_r{i}(j,1:max(length(find(m_loc{j-tFIR}(k,:))))));
            end
        end
    else
        for j = 1:tFIR
            Psi_rows{i} = Psi(r{i},s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))));
            Lambda_rows{i} = Lambda(r{i},s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))));
        end
    end
end        

end