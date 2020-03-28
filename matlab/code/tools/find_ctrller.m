function ctrller = find_ctrller(sys, clmaps, params, eqnErrCoeff)
% Find alternate implementation, returned in ctrller
% Returns
%    ctrller     : Ctrller containing implementation matrices (Rc, Mc)
% Inputs
%    sys         : LTISystem containing system matrices
%    clmaps      : contains closed loop maps (R, M)
%    params      : SLSParams containing parameters
%    eqnErrCoeff : regularization coefficient on equation error term
fprintf('Finding controller\n\n');
params.print()
params.approx_ = true;
params.sanity_check()
    
Tc = params.T_;

cvx_begin % TODO: quiet?

expression Rcs(sys.Nx, sys.Nx, Tc)
expression Mcs(sys.Nu, sys.Nx, Tc)
expression Deltas(sys.Nx, sys.Nx, Tc)

% populate decision variables for ease-of-use
Rc    = cell(Tc, 1); 
Mc    = cell(Tc, 1);
Delta = cell(Tc, 1);
for t = 1:Tc
    Rc{t} = Rcs(:,:,t); Mc{t} = Mcs(:,:,t); Delta{t} = Deltas(:,:,t);
end

% enforce constraints
if ~isempty(params.constraints_)
    [RSupp, MSupp, count] = get_supports(sys, params);
    variable RMcSuppVals(count)
    [Rc, Mc] = add_sparse_constraints(Rc, Mc, RSupp, MSupp, RMcSuppVals, Tc);
end

objective = get_total_objective(sys, params, Rc, Mc);
objective = objective + get_eqn_err_obj(sys, clMaps, Rc, Mc);

Rc{1} == eye(sys.Nx);

for t=1:Tc-1
    Delta{t} = Rc{t+1} - sys.A*Rc{t} - sys.B2*Mc{t};
end
Delta{Tc} = - sys.A*Rc{Tc} - sys.B2*Mc{Tc};
% regularization for stability
objective = objective + params.approxCoeff_ * get_stab_obj(Delta);

minimize(objective);
cvx_end

% outputs
ctrller              = Ctrller;
ctrller.Rc_          = Rc;
ctrller.Mc_          = Mc;
ctrller.solveStatus_ = cvx_status;

if strcmp(cvx_status, 'Solved')
    fprintf(['Solved!', '\n\n']);
else
    sls_warning(['Solver exited with status ', cvx_status]);
end
end
