function [R, M] = add_sparse_constraints(R, M, RSupp, MSupp, X, T)
    spot = 0;
    for t = 1:T
        suppSizeR = sum(sum(RSupp{t}));
        R{t}(RSupp{t}) = X(spot+1:spot+suppSizeR);
        spot = spot + suppSizeR;

        suppSizeM = sum(sum(MSupp{t}));
        M{t}(MSupp{t}) = X(spot+1:spot+suppSizeM);
        spot = spot + suppSizeM;
    end
end