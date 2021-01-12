function r = assign_rows(sys, tFIR)
% r{i}   represents the set of rows subsystem i solves for

Nx = sys.Nx; Nu = sys.Nu;
r  = cell(Nx, 1);

for i = 1:Nx
    actuator = find(sys.B2(i,:)); % assumes 1 actuator per system
    
    nRows = tFIR;
    if ~isempty(actuator) % this subsystem has actuation
        nRows    = nRows + tFIR-1; % include rows from Phi_u
    end
    
    r{i} = cell(nRows, 1);
    
    for j=1:tFIR 
        r{i}{j}   = (j-1)*Nx + i; % rows representing state i
    end
    if ~isempty(actuator)
        for j=1:tFIR-1
            j_ = j + tFIR;
            r{i}{j_}   = tFIR*Nx + (j-1)*Nu + actuator;
        end
    end
end

end