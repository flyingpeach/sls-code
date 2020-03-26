function total_objective = get_total_objective(sys, params, R, M)

total_objective = 0;

for i=1:length(params.objectives_)
    this_obj  = params.objectives_{i};
    objType   = this_obj{1};
    objRegVal = this_obj{2};

    objective = 0;
    switch objType
        case SLSObjective.H2
            objective = get_H2_obj(sys, R, M);
        case SLSObjective.HInf
            objective = get_HInf_obj(sys, R, M);
        case SLSObjective.L1
            objective = get_L1_obj(sys, R, M);
        case SLSObjective.RFD
            objective = get_rfd_obj(sys, M);
    end
    
    total_objective = total_objective + objRegVal * objective;    
end

end