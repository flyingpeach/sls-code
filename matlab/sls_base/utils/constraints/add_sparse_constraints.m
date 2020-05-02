function [R, M] = add_sparse_constraints(R, M, RSupp, MSupp, X, T)
    spot = 0;
    for t = 1:T
        suppR = find(RSupp{t});
        num = sum(sum(RSupp{t}));
        R{t}(suppR) = X(spot+1:spot+num);
        spot = spot + num;

        suppM = find(MSupp{t});
        num = sum(sum(MSupp{t}));
        M{t}(suppM) = X(spot+1:spot+num);
        spot = spot + num;
    end
end