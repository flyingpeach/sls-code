function r = assign_rows_h_terminal(sys, params)
% r{i} represents the set of rows subsystem i solves for
% When terminal set constraints are present, use this function instead
% of assign_rows_h
% Also contains sanity check for locality for terminal set 

r = cell(Nx, 1);
% Assign rows of original H matrix (minus state constraint at time T)
% assuming UB/LB on all states, and no input bounds
for i=1:Nx
    for t=1:T-1
        r{i}(end+1) = Nx*(t-1) + i;
    end
end

HT         = params.terminal_H_;
Comms_Adj  = abs(sys.A) > 0;
StateSupp  = Comms_Adj^(params.locality_-1) > 0;

for rowHT=1:size(HT, 1)    
    Comms_Hi    = HT(rowHT, :) ~= 0; % boolean
    locality_ok = false;
    rowH        = rowHT + 2*Nx*(T-1);
    
    for i=find(HT(rowHT, :))
        if isempty(find(StateSupp(i, :) - Comms_Hi < 0, 1))
            % H(i,:) obeys locality sparsity for some node
            % We can say H(i,:) "belongs" to this node
            r{i}(end+1) = rowH;
            locality_ok = true;
        end
    end
    
    if locality_ok == false
        mpc_error('Terminal set coupling violates locality!');
    end
end

end