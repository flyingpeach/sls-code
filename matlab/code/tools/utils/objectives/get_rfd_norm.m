function objective = get_rfd_norm(sys, M)
objective = 0;
for i = 1:sys.Nu
    Mi = [];
    for t = 1:length(M)
        Mi = [Mi, M{t}(i,:)];
    end    
    objective = objective + norm(Mi, 2);
end
end