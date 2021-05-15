function r = assign_rows_phi(sys, tFIR)
% r{i}   represents the set of rows subsystem i solves for

Nx = sys.Nx; Nu = sys.Nu;
r  = cell(Nx, 1);

for i = 1:Nx
    actuator = find(sys.B2(i,:)); % assumes 1 actuator per system
    r{i} = [];
    
    for j=1:tFIR 
        r{i}(end+1) = (j-1)*Nx + i; % rows representing state i
    end
    if ~isempty(actuator)
        for j=1:tFIR-1
            r{i}(end+1) = tFIR*Nx + (j-1)*Nu + actuator;
        end
    end
end

end