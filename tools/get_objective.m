function objective = get_objective(sys, params, R, M)
% Gets the objective function value based on the 
% Returns 
%    objective : objective function value
% Inputs
%    sys       : an LTISystem
%    params    : SLSParams containing parameters

switch params.obj_
    case Objective.TrajTrack
        objective = traj_track(sys, params, R, M);
    case Objective.H2
        objective = compute_H2(sys, params, R, M);
    case Objective.HInf
        objective = compute_Hinf(sys, params, R, M);
    otherwise
        objective = 0;
        disp('[SLS WARNING] Objective = constant, only finding feasible solution')
end
end


% local functions
function objective = traj_track(sys, params, R, M)
% Reformulation of (x-xd)' Q (x-xd) + u' R u
% i.e. track desired trajectory while minimizing input

% TODO: should be able to pass these in
xPenalty = eye(sys.Nx); % matrix Q in cost function
uPenalty = 0.1 * eye(sys.Nu); % matrix R in cost function

% TODO: temp hack for proof-of-concept; will need to max over w
w = ones(sys.Nx, 40);

objective = 0;
for t = 1:params.tFIR_
    % need to do the vect operation because of quirk in cvx
    vect = vec(blkdiag(sqrtm(xPenalty), sqrtm(uPenalty)) * ...
           ([R{t}; M{t}] * w(:,t) - [params.xDes_(:,t); zeros(sys.Nu, 1)]));

    objective = objective + vect'*vect;
end
end


function objective = compute_H2(sys, params, R, M)
% return ||[C1,D12][R;M]||_H2^2 as per (4.20)

objective = 0;
for t = 1:params.tFIR_
    %need to do the vect operation because of quirk in cvx
    vect = vec([sys.C1, sys.D12]*[R{t};M{t}]);
    objective = objective + vect'*vect;
end
end


function objective = compute_Hinf(sys, params, R, M)
% return max singular value of [C1,D12][R;M]

mtx = [];
for t = 1:params.tFIR_
    mtx = blkdiag(mtx, [sys.C1, sys.D12]*[R{t};M{t}]);
end

objective = sigma_max(mtx);
end