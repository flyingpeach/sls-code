function [x, u, avgTime, avgIter] = sls_mpc(sys, x0, params)

    avgIter = [];
    
    if params.mode_ == MPCMode.Centralized
        [x, u, avgTime] = mpc_centralized(sys, x0, params);
    elseif params.mode_ == MPCMode.Distributed
        % Implementation of distributed MPC currently requires
        %  1. One actuator per subsystem
        %  2. Actuator indices should be a strictly increasing
        %     function of the subsystem indices they correspond to
        %  Example: B2 = [0 0 1;
        %                [0 1 0] violates (2)
        %           since actuator 3 corresponds to subsystem 1
        %           but   actuator 2 corresponds to subsystem 2
        %  B2 = [0 1 0;
        %       [0 0 1] will work
        
        lastActIdx = 0;
        for row=1:sys.Nu
            if length(find(sys.B2(row,:))) > 1
                mpc_error('Maximum one actuator per subsystem');
            end
            if find(sys.B2(row,:)) <= lastActIdx
                mpc_error('sys.B2 has unsupported sparsity pattern');
            end
            lastActIdx = find(sys.B2(row,:));
        end
     
        if ~params.has_coupling()
            [x, u, avgTime, avgIter] = mpc_algorithm_1(sys, x0, params);
        else
            % TODO: check cost/constraint matrices should be locality-limited      
            [x, u, avgTime, avgIter] = mpc_algorithm_2(sys, x0, params);
        end
    else
        mpc_error('Unrecognized MPC mode specified!');
    end    
    
end