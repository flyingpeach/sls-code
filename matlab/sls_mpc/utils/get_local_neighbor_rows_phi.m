function r = get_local_neighbor_rows_phi(sys, params, rPhi)
% r: per processor, the rows that belong to other processors that this
% processor can access (e.g. other processors within this local patch)

% constraint matrix 
Nx = sys.Nx;
r  = cell(Nx, 1);

Comms_Adj = abs(sys.A)>0;
RSupp     = Comms_Adj^(params.locality_-1) > 0;

for i = 1:Nx
    neighbors = find(RSupp(i, :));
    for neighbor = neighbors
        if neighbor ~= i % don't count self
            
            % we have access to all neighbors' rows
            r{i} = [r{i} rPhi{neighbor}];            
        end
    end
end

end
