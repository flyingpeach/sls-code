function [r, s_r] = get_row_locality(sys, tFIR, r_loc, m_loc)

Nx = sys.Nx; Nu = sys.Nu;

r     = cell(1, Nx);
s_r   = cell(1, Nx);

k = 0;
for i = 1:Nx
    if mod(i, Nx/Nu) == 0 % Decide whether or not there is actuation
        s_r{i} = zeros(tFIR+(tFIR-1),Nx); % Preallocate the indices
        k = k+1;
        for j = 1:tFIR+(tFIR-1)
            if j<=tFIR
                r{i}(j) = Nx*(j-1) + i;
                s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))) = find(r_loc{j}(i,:));
            else
                r{i}(j) = Nu*(j-tFIR-1) + Nx*tFIR + k;
                s_r{i}(j,1:max(length(find(m_loc{j-tFIR}(k,:))))) = find(m_loc{j-tFIR}(k,:));
            end
        end
    else
        s_r{i} = zeros(tFIR,Nx); % Preallocate the indices
        for j = 1:tFIR
            r{i}(j) = Nx*(j-1) + i;
            s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))) = find(r_loc{j}(i,:));
        end
    end
    s_r{i}( :, ~any(s_r{i},1) ) = []; % Eliminate columns with only zeros
end