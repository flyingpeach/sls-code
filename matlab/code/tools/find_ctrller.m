function ctrller = find_ctrller(sys, clMaps, params, eqnErrCoeff)
% Find alternate implementation, returned in ctrller
% Returns
%    ctrller     : Ctrller containing implementation matrices (Rc, Mc)
% Inputs
%    sys         : LTISystem containing system matrices
%    clMaps      : contains closed loop maps (R, M)
%    params      : SLSParams containing parameters
%    eqnErrCoeff : regularization coefficient on equation error term
fprintf('Finding controller\n\n');
params.approx_ = true;
params.add_objective(SLSObjective.EqnErr, eqnErrCoeff);

slsOut  = state_fdbk_sls(sys, params, clMaps);
ctrller = Ctrller.ctrller_from_cl_maps(slsOut);
end
