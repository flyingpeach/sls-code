function [xs, us, avgStats] = sls_mpc(sys, x0, params, tHorizon)
% Call this function once to run mpc on a non-time-varying system
% Wrapper function for mpc_distributed and mpc_centralized, which are 
% per-timestep functions
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu 
%   x0      : Starting system state
%   params  : MPCParams containing parameters for mpc
%   tHorizon: total time to run MPC for
% Outputs
%   avgStats: average runtime, iters, consIters (per inner loop iter) per timestep
%             iters and consIters will be empty for centralized MPC

xs      = zeros(sys.Nx, tHorizon);
us      = zeros(sys.Nu, tHorizon);
xs(:,1) = x0;

times     = zeros(tHorizon-1, 1);
iters     = zeros(tHorizon-1, 1);
consIters = zeros(tHorizon-1, 1);

% initially no warm start; will be populated for subsequent timesteps
warmStart = [];

for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t+1, tHorizon);
    
    if params.mode_ == MPCMode.Distributed
        if params.accounts_for_disturbance() % robust MPC
            [xs(:,t+1), us(:,t), stats, warmStart] = rmpc_distributed(sys, xs(:,t), params, warmStart);
        else % noiseless MPC
            [xs(:,t+1), us(:,t), stats, warmStart] = mpc_distributed(sys, xs(:,t), params, warmStart);
        end        
        times(t)     = stats.time_;
        iters(t)     = stats.iters_;
        consIters(t) = stats.consIters_;
        
    elseif params.mode_ == MPCMode.Centralized % centralized algorithms return no iteration info
        if params.accounts_for_disturbance() % robust MPC
            [xs(:,t+1), us(:,t), times(t)] = rmpc_centralized(sys, xs(:,t), params);
        else % standard MPC
            [xs(:,t+1), us(:,t), times(t)] = mpc_centralized(sys, xs(:,t), params);
        end
        
    else
        mpc_error('Unrecognized MPC mode specified!');
    end
            
end

avgStats       = MPCStats();
avgStats.time_ = mean(times);
if params.mode_ == MPCMode.Distributed
    avgStats.iters_     = mean(iters);
    avgStats.consIters_ = mean(consIters);
end

end