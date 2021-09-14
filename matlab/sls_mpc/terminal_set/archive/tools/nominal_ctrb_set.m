function [H_, h_] = nominal_ctrb_set(sys, params, N)
% Gives N-step controllable set for some set defined by
% state constraints in params. params also contains input constraints
% Output set is defined as x: H_*x <= h_

T = N+1;

tFIR         = params.tFIR_;
params.tFIR_ = T; % temporarily change param (for get_constraint functions)

ZAB        = get_sls_constraint(sys, T);
[Hxu, hxu] = get_sys_constraints(sys, params);
[H_, h_]   = nominal_set(Hxu, ZAB, hxu, sys.Nx);

% restore parameter
params.tFIR_ = tFIR;

end