function clMaps = state_fdbk_sls(sys, params)
% System level synthesis with state feedback
% Returns 
%    slsOuts: SLSOutputs containing system responses and other info
% Inputs
%    sys     : LTISystem containing system matrices
%    params  : SLSParams containing parameters
params.print()
params.sanity_check()

cvx_begin quiet

if isempty(params.constraints_) % SLS with no constraints
    variable Rs(sys.Nx, sys.Nx, params.T_)
    variable Ms(sys.Nu, sys.Nx, params.T_)
else 
    expression Rs(sys.Nx, sys.Nx, params.T_)
    expression Ms(sys.Nu, sys.Nx, params.T_)
end

if params.approx_
    expression Delta(sys.Nx, sys.Nx * params.T_)
end

% populate decision variables for ease-of-use
R = cell(params.T_, 1); 
M = cell(params.T_, 1);
for t = 1:params.T_
    R{t} = Rs(:,:,t); M{t} = Ms(:,:,t);
end

% enforce constraints
if ~isempty(params.constraints_)
    [R, M] = add_sparse_constraints(sys, params, R, M);
end
objective = get_total_objective(sys, params, R, M);

% achievability  / approx achievability constraints
R{1} == eye(sys.Nx);
R{params.T_} == zeros(sys.Nx, sys.Nx);

if params.approx_
    for t=1:params.T_-1
        Delta(:,(t-1)*sys.Nx+1:t*sys.Nx) = R{t+1} - sys.A*R{t} - sys.B2*M{t};
    end
    objective = objective + params.approxCoeff_ * norm(Delta, inf);    
else
    for t=1:params.T_-1
        R{t+1} == sys.A*R{t} + sys.B2*M{t};
    end
end

minimize(objective);
cvx_end

% outputs
slsOuts              = CLMaps();
slsOuts.acts_        = get_acts_rfd(sys, params, M); % rfd actuator selection
slsOuts.R_           = R;
slsOuts.M_           = M;
slsOuts.solveStatus_ = cvx_status;

if strcmp(cvx_status, 'Solved')
    fprintf(['Solved!', '\n']);
else
    sls_warning(['Solver exited with status ', cvx_status]);
end
end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acts = get_acts_rfd(sys, params, M)
tol = 1e-4;

acts = [];
if params.rfd_
    for i=1:sys.Nu
        if norm(vec(M{1}(i,:)),2) > tol
            acts = [acts; i];
        end
    end    
end
end