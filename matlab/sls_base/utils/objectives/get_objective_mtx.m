function objectiveMtx = get_objective_mtx(sys, R, M)

T = length(R);
objectiveMtx = cell(T, 1);
for k=1:T
    objectiveMtx{k} = [sys.C1, sys.D12]*[R{k};M{k}];
end

end