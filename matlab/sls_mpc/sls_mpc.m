function [xs, us, avgTime, avgIters] = sls_mpc(sys, x0, params, tHorizon)
% Call this function once to run mpc on a system with 
% non time-varying parameters / dynamics
% Wrapper function for sls_mpc_per_time
% Allows some time-varying systems and constraints
%   avgTime  : Runtime per state per timestep
%   avgIters : Total ADMM iters per state per timestep (for distributed MPC only)
% For time / iteration calculations, t=1 is omitted to omit warm-up effects

xs      = zeros(sys.Nx, tHorizon);
us      = zeros(sys.Nu, tHorizon);
xs(:,1) = x0;

times = zeros(tHorizon-1, 1);
iters = zeros(tHorizon-1, 1);

for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t, tHorizon-1);
    if params.mode_ == MPCMode.Distributed
        [xs(:,t+1), us(:,t), times(t), iters(t)] = sls_mpc_per_time(sys, xs(:,t), params);
    else % no iters in centralized MPC
        [xs(:,t+1), us(:,t), times(t)] = sls_mpc_per_time(sys, xs(:,t), params);
    end
end

avgTime = mean(times(2:end)); % omit t=1

if params.mode_ == MPCMode.Distributed
    avgIters = mean(iters(2:end));
else
    avgIters = [];
end

end