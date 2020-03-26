function total_objective = get_total_objective(sys, params, R, M)

total_objective = 0;

for i=1:length(params.objectives_)
    this_obj  = params.objectives_{i};
    objType   = this_obj(1);
    objRegVal = this_obj(2);

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



% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function objective = get_H2_obj(sys, R, M)
objective = 0;
for t = 1:length(R)
    % need to do the vect operation because of quirk in cvx
    vect = vec([sys.C1, sys.D12]*[R{t};M{t}]);
    objective = objective + vect'*vect;
end
end


function objective = get_HInf_obj(sys, R, M)
mtx = [];
for t = 1:length(R)
    mtx = blkdiag(mtx, [sys.C1, sys.D12]*[R{t};M{t}]);
end

objective = sigma_max(full(mtx));
end


function objective = get_L1_obj(sys, R, M)
mtx = [];
for t = 1:length(R)
    mtx = blkdiag(mtx, [sys.C1, sys.D12]*[R{t};M{t}]);
end

objective = norm(mtx, Inf); % note: L1 is induced inf-inf norm
end


function objective = get_rfd_obj(sys, M)
objective = 0;
for i = 1:sys.Nu
    Mi = [];
    for t = 1:length(M)
        Mi = [Mi, M{t}(i,:)];
    end    
    objective = objective + norm(Mi, 2);
end
end