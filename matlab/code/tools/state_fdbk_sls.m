function clMaps = state_fdbk_sls(sys, params)
% System level synthesis with state feedback
% Returns 
%    clMaps : SLSOutputs containing system responses and other info
% Inputs
%    sys     : LTISystem containing system matrices
%    params  : SLSParams containing parameters
fprintf('Finding closed loop maps\n\n');
params.print()
params.sanity_check()

T = params.T_;

cvx_begin quiet

if isempty(params.constraints_) % SLS with no constraints
    variable Rs(sys.Nx, sys.Nx, T)
    variable Ms(sys.Nu, sys.Nx, T)
else 
    expression Rs(sys.Nx, sys.Nx, T)
    expression Ms(sys.Nu, sys.Nx, T)
end

if params.approx_
    expression Delta(sys.Nx, sys.Nx * T)
end

% populate decision variables for ease-of-use
R = cell(T, 1); 
M = cell(T, 1);
for t = 1:T
    R{t} = Rs(:,:,t); M{t} = Ms(:,:,t);
end

% enforce constraints
if ~isempty(params.constraints_)
    [RSupp, MSupp, count] = get_supports(sys, params);
    variable RMSuppVals(count)
    [R, M] = add_sparse_constraints(R, M, RSupp, MSupp, RMSuppVals, T);
end
objective = get_total_objective(sys, params, R, M);

% achievability  / approx achievability constraints
R{1} == eye(sys.Nx);
R{T} == zeros(sys.Nx, sys.Nx);

if params.approx_
    for t=1:T-1
        Delta(:,(t-1)*sys.Nx+1:t*sys.Nx) = R{t+1} - sys.A*R{t} - sys.B2*M{t};
    end
    objective = objective + params.approxCoeff_ * norm(Delta, inf);    
else
    for t=1:T-1
        R{t+1} == sys.A*R{t} + sys.B2*M{t};
    end
end

minimize(objective);
cvx_end

% outputs
clMaps              = CLMaps();
clMaps.acts_        = get_acts_rfd(sys, M); % rfd actuator selection
clMaps.R_           = R;
clMaps.M_           = M;
clMaps.solveStatus_ = cvx_status;

if strcmp(cvx_status, 'Solved')
    fprintf(['Solved!', '\n\n']);
else
    sls_warning(['Solver exited with status ', cvx_status]);
end
end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acts = get_acts_rfd(sys, M)
tol = 1e-4;

acts = [];
    for i=1:sys.Nu
        if norm(vec(M{1}(i,:)),2) > tol
            acts = [acts; i];
        end
    end
end