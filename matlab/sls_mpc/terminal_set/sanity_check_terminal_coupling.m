function sanity_check_terminal_coupling(sys, params)
% Note: This only checks state coupling, not input coupling

H = params.terminal_H_;

Comms_Adj  = abs(sys.A) > 0;
StateSupp  = Comms_Adj^(params.locality_-1) > 0;

for i=1:size(H, 1)
    
    Comms_Hi = H(i, :) ~= 0; % boolean
    nodes = find(H(i, :));
    
    locality_ok = false;
    for node=nodes        
        if isempty(find(StateSupp(node, :) - Comms_Hi < 0, 1))
            % H(i,:) obeys locality sparsity for some node
            % We can say H(i,:) "belongs" to this node
            locality_ok = true;
            break
        end
    end
    
    if locality_ok == false
        mpc_error('Terminal set coupling violates locality!');
    end
end

end