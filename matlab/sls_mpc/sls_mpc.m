function [x, u, avgTime, avgIter] = sls_mpc(sys, x0, params)

    if ~params.has_coupling()
        [x, u, avgTime, avgIter] = mpc_algorithm_1(sys, x0, params);
    else
        [x, u, avgTime, avgIter] = mpc_algorithm_2(sys, x0, params);
    end

end