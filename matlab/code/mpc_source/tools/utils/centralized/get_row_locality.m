function [r, s_r] = get_row_locality(sys, tFIR, r_loc, m_loc)

Nx = sys.Nx; Nu = sys.Nu;

r     = cell(1, Nx);
s_r   = cell(1, Nx);

% TODO: this only works for quasi-diagonal B2 matrix
actuator = 0;

for i = 1:Nx
    hasActuation = ~all(sys.B2(i,:) == 0);
    
    nRows = tFIR;
    if hasActuation
        nRows    = nRows + tFIR-1; % include rows from M
        actuator = actuator + 1;
    end
    
    r{i}   = zeros(1, nRows);
    s_r{i} = cell(nRows, 1);
    
    for j=1:tFIR
        r{i}(j)   = (j-1)*Nx + i;       
        s_r{i}{j} = find(r_loc(i, :));
    end
    if hasActuation
        for j=1:tFIR-1
            j_ = j + tFIR;
            r{i}(j_)   = tFIR*Nx + (j-1)*Nu + actuator;
            s_r{i}{j_} = find(m_loc(actuator, :));
        end
    end
end

end