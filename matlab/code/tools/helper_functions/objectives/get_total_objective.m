function totalObjective = get_total_objective(sys, params, R, M)

totalObjective = 0;

objectiveMtx = cell(params.T_, 1);
for k=1:params.T_
    objectiveMtx{k} = [sys.C1, sys.D12]*[R{k};M{k}];
end

for i=1:length(params.objectives_)
    this_obj  = params.objectives_{i};
    objType   = this_obj{1};
    objRegVal = this_obj{2};

    objective = 0;
    switch objType
        case SLSObjective.H2
            objective = get_H2_obj(objectiveMtx);
        case SLSObjective.HInf
            objective = get_HInf_obj(objectiveMtx);
        case SLSObjective.L1
            objective = get_L1_obj(objectiveMtx);
        case SLSObjective.RFD
            objective = get_rfd_obj(sys, M);
    end
    
    totalObjective = totalObjective + objRegVal * objective;    
end

end