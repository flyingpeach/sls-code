function slsOuts = state_fdbk_sls(sys, params)
% System level synthesis with state feedback
% Can optionally use regularization for design (rfd) as well
% Returns 
%    slsouts: SLSOutputs containing system responses and other info
% Inputs
%    sys    : LTISystem containing system matrices
%    params : SLSParams containing parameters

slsOuts = SLSOutputs;

% used in rfd
tol = 1e-4;

if params.mode_ ~= SLSMode.Basic  
    make_d_localized_constraints(sys, params);
    [Rsupport, Msupport, count] = make_d_localized_constraints(sys, params);
end

cvx_begin

% decision variables
if params.mode_ ~= SLSMode.Basic
    variable X(count)
    expression Rs(sys.Nx, sys.Nx, params.tFIR_)
    expression Ms(sys.Nu, sys.Nx, params.tFIR_)
    if params.mode_ == SLSMode.ApproxDLocalized
        expression Delta(sys.Nx, sys.Nx * params.tFIR_)
    end
else % basic SLS
    variable Rs(sys.Nx, sys.Nx, params.tFIR_)
    variable Ms(sys.Nu, sys.Nx, params.tFIR_)
end

% populate decision variables
for t = 1:params.tFIR_
    R{t} = Rs(:,:,t);
    M{t} = Ms(:,:,t);
end

% locality constraints 
% automatically enforced by limiting support of R, M
if params.mode_ ~= SLSMode.Basic  
    spot = 0;
    for t = 1:params.tFIR_
        suppR = find(Rsupport{t});
        num = sum(sum(Rsupport{t}));
        R{t}(suppR) = X(spot+1:spot+num);
        spot = spot + num;

        suppM = find(Msupport{t});
        num = sum(sum(Msupport{t}));
        M{t}(suppM) = X(spot+1:spot+num);
        spot = spot + num;
    end
end

objective   = get_objective(sys, params, R, M);
robust_stab = 0;

% rfd actuation penalty
act_penalty = 0;
if params.rfd_
    for i = 1:sys.Nu
        Mi = [];
        for t = 1:params.tFIR_
            Mi = [Mi, M{t}(i,:)];
        end    
        act_penalty = act_penalty + norm(Mi, 2);
    end
end

% achievability  / approx achievability constraints
R{1} == eye(sys.Nx);
R{params.tFIR_} == zeros(sys.Nx, sys.Nx);

if params.mode_ == SLSMode.ApproxDLocalized
    for t=1:params.tFIR_-1
        Delta(:,(t-1)*sys.Nx+1:t*sys.Nx) = R{t+1} - sys.A*R{t} - sys.B2*M{t};
    end
    robust_stab = norm(Delta, inf); % < 1 means we can guarantee stab
    objective = objective + params.robCoeff_ * robust_stab;    
else
    for t=1:params.tFIR_-1
        R{t+1} == sys.A*R{t} + sys.B2*M{t};
    end
end

if params.rfd_
    objective = objective + params.rfdCoeff_ * act_penalty;
end

% solve minimization problem
minimize(objective);
cvx_end

% rfd actuator selection
if params.rfd_
    acts = [];
    for i=1:sys.Nu
        if norm(vec(M{1}(i,:)),2) > tol
            acts = [acts; i];
        end
    end
    slsOuts.acts_ = acts;
else
    % not set in case of rfd
    % doesn't include regularization terms
    slsOuts.clnorm_ = get_objective(sys, params, R, M);
end

slsOuts.R_           = R;
slsOuts.M_           = M;
slsOuts.robust_stab_ = robust_stab;