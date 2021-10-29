function r = assign_rows_h_terminal_only(sys, params, HT)
% This only assign rows for the terminal state constraint matrix
% Also contains sanity check for locality for terminal set 

r = cell(sys.Nx, 1);

commsAdj  = abs(sys.A) > 0;
stateSupp = commsAdj^(params.locality_-1) > 0;

for rowHT=1:size(HT, 1)    
    Comms_Hi    = HT(rowHT, :) ~= 0; % boolean
    locality_ok = false;
    
    if isempty(find(HT(rowHT,:)))
        % Empty row; add to first subsystem by default
        r{1}(end+1) = rowHT; 
        locality_ok = true;
    end
    
    for i=find(HT(rowHT, :))
        if isempty(find(stateSupp(i, :) - Comms_Hi < 0, 1))
            % H(i,:) obeys locality sparsity for some node
            % We can say H(i,:) "belongs" to this node
            r{i}(end+1) = rowHT;
            locality_ok = true;
            break;
        end
    end
    
    if locality_ok == false
        mpc_error('Terminal set coupling violates locality!');
    end
end

end