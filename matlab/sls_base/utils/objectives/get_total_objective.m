function totalObjective = get_total_objective(sys, params, R, M, varargin)
clMapsIn = [];
try 
    clMapsIn = varargin{1};
end

totalObjective = 0;
objectiveMtx   = get_objective_mtx(sys, R, M);

for i=1:length(params.objectives_)
    this_obj  = params.objectives_{i};
    objType   = this_obj{1};
    objRegVal = this_obj{2};

    objective = 0;
    switch objType
        case SLSObjective.H2
            objective = get_H2_norm(objectiveMtx);
        case SLSObjective.HInf
            objective = get_HInf_norm(objectiveMtx);
        case SLSObjective.L1
            objective = get_L1_norm(objectiveMtx);
        case SLSObjective.Eps1
            objective = get_eps1_norm(objectiveMtx);
        case SLSObjective.OneToOne
            objective = get_1to1_norm(objectiveMtx);
        case SLSObjective.RFD
            objective = get_rfd_norm(sys, M);
        case SLSObjective.EqnErr
            objective = get_eqn_err(sys, clMapsIn, R, M);
    end
    
    totalObjective = totalObjective + objRegVal * objective;    
end

end