function [x, u, avgTime, avgIter] = sls_mpc(sys, x0, params)

    avgIter = [];
    
    if params.mode_ == MPCMode.Centralized
        [x, u, avgTime] = mpc_centralized(sys, x0, params);
    elseif params.mode_ == MPCMode.Distributed
        if ~params.has_coupling()
            [x, u, avgTime, avgIter] = mpc_algorithm_1(sys, x0, params);
        else
            [x, u, avgTime, avgIter] = mpc_algorithm_2(sys, x0, params);
        end
    else
        mpc_error('Unrecognized MPC mode specified!');
    end    
    
end