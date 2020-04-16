function [c, s_c] = get_column_locality(sys, tFIR, r_loc, m_loc)

Nx = sys.Nx; Nu = sys.Nu;

c     = cell(1, Nx);
s_c   = cell(1, Nx);

for i = 1:Nx
    c{i} = i;
    count = 0;
    for j = 1:tFIR+(tFIR-1)
        if j<=tFIR
            find_locR = find(r_loc{j}(:,i));
            for k =1:max(length(find_locR))
                count = count +1;
                s_c{i}(count) = find_locR(k)+(j-1)*Nx;
                if j == tFIR
                end
            end
        else
            find_locM = find(m_loc{j-tFIR}(:,i));
            for k =1:max(length(find_locM))
                count = count +1;
                s_c{i}(count) = find_locM(k)+(j-tFIR-1)*Nu+tFIR*Nx;
            end
        end
    end
end