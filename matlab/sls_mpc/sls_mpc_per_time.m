function [x, u, time, iters] = sls_mpc_per_time(sys, x0, params)
% Call this function once per time-step with the desired system / params
% Allows some time-varying systems and constraints
%
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu 
%   x0      : Starting system state
%   params  : MPCParams containing parameters for mpc
% Outputs
%   x       : Next state if MPC input is used
%   u       : Current input as calculated by MPC
%   time    : Total runtime per state
%   iters   : Total ADMM iters per state (for distributed MPC only)

iters = [];

if params.mode_ == MPCMode.Centralized
    if params.accounts_for_disturbance() % robust MPC    
        [x, u, time] = rmpc_centralized(sys, x0, params);
    else
        [x, u, time] = mpc_centralized(sys, x0, params);
    end
elseif params.mode_ == MPCMode.Distributed
    if params.accounts_for_disturbance() % robust MPC
        [x, u, time, iters] = rmpc_distributed(sys, x0, params);
    else
        [x, u, time, iters] = mpc_distributed(sys, x0, params);
    end        
else
    mpc_error('Unrecognized MPC mode specified!');
end    
end