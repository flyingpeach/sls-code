function [xs, us, avgTime, avgIters, avgConsIters] = sls_mpc(sys, x0, params, tHorizon)
% Call this function once to run mpc on a non-time-varying system
% Wrapper function for mpc_distributed and mpc_centralized, which are 
% per-timestep functions
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu 
%   x0      : Starting system state
%   params  : MPCParams containing parameters for mpc
%   tHorizon: total time to run MPC for
% Outputs
%   avgTime     : Runtime per state per timestep
%   avgIters    : ADMM iters per state per timestep (distributed MPC only)
%   avgConsIters: ADMM consensus iters per state, per inner loop iter, per timestep (distributed MPC only)
% For time / iteration calculations, t=1 is omitted to omit warm-up effects

xs      = zeros(sys.Nx, tHorizon);
us      = zeros(sys.Nu, tHorizon);
xs(:,1) = x0;

times     = zeros(tHorizon-1, 1);
iters     = zeros(tHorizon-1, 1);
consIters = zeros(tHorizon-1, 1);

for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t+1, tHorizon);
    if params.mode_ == MPCMode.Distributed
        if params.accounts_for_disturbance() % robust MPC
            [xs(:,t+1), us(:,t), times(t), iters(t)] = rmpc_distributed(sys, xs(:,t), params);
        else % standard MPC
            [xs(:,t+1), us(:,t), times(t), iters(t), consIters(t)] = mpc_distributed(sys, xs(:,t), params);
        end
        
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

avgTime = mean(times(2:end)); % omit t=1

avgIters     = [];
avgConsIters = [];
if params.mode_ == MPCMode.Distributed
    avgIters = mean(iters(2:end));
    if ~params.accounts_for_disturbance
        avgConsIters = mean(consIters(2:end));
    end
end

end