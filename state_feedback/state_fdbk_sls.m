function [R, M, clnorm, robust_stab] = state_fdbk_sls(sys, params)
% System level synthesis with state feedback

% The approx d-localized version currently only enforces l1->l1 small gain 
% on Delta; TODO: will need to generalize

% Relies on robust variant of theory introduced in 
% http://www.its.caltech.edu/~james/papers/v_localizable.pdf

% Returns 
%    R, M   : system response as defined in Thm 1 of 
%             https://arxiv.org/pdf/1610.04815.pdf
%    clnorm : final (optimal value) of objective
% Inputs
%    sys    : an LTISystem
%    params : SLSParams containing parameters

if params.mode_ ~= SLSMode.Basic  
    make_d_localized_constraints(sys, params);
    [Rsupport, Msupport, count] = make_d_localized_constraints(sys, params);
end

cvx_begin
cvx_precision low

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

objective = get_objective(sys, params, R, M);

% achievability  / approx achievability constraints
R{1} == eye(sys.Nx);

if params.mode_ == SLSMode.ApproxDLocalized
    for t=1:params.tFIR_-1
        Delta(:,(t-1)*sys.Nx+1:t*sys.Nx) = R{t+1} - sys.A*R{t} - sys.B2*M{t};
    end
    robust_stab = norm(Delta, inf); % < 1 means we can guarantee stab
    minimize(objective + params.lambda_ * robust_stab);

else
    for t=1:params.tFIR_-1
        R{t+1} == sys.A*R{t} + sys.B2*M{t};
    end
    R{params.tFIR_} == zeros(sys.Nx, sys.Nx);
    minimize(objective);
end

cvx_end
clnorm = objective;