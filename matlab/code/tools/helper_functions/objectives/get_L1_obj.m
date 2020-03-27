function objective = get_L1_obj(sys, R, M)
objective = 0;
for t = 1:length(R)
    objective = max(objective, norm([sys.C1, sys.D12]*[R{t};M{t}], Inf));
end
end