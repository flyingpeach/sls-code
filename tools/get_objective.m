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
        warning('Objective = constant, only finding feasible solution')
end