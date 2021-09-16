function [cs, ms, gs] = precalculate_col_robust(sys, T, Cost, H, s_cPsi, s_cLambda)

Nx   = sys.Nx;
nPhi = Nx*T + sys.Nu*(T-1);

ZAB = get_constraint_zab(sys, T);
Eye = eye(Nx*T);

cs = cell(Nx*T, 1);
ms = cell(Nx*T, 1);
gs = cell(Nx*T, 1);

M1 = [Cost; H];
for col = 1:length(s_cPsi)        
    zab_     = ZAB(:, s_cPsi{col});
    zeroRows = find(all(zab_ == 0, 2));
    keepRows = setdiff(1:Nx*T, zeroRows);
    f_       = ZAB(keepRows, s_cPsi{col});
    gs{col}  = Eye(keepRows, col);
        
    if col <= Nx
        m_ = M1(s_cLambda{col}, s_cPsi{col});
    else
        m_ = H(s_cLambda{col}-nPhi, s_cPsi{col});
    end
    n = size(m_, 2);
    p = size(f_, 1);

    c       = pinv([m_'*m_ f_'; f_ zeros(p,p)]);
    cs{col} = c(1:n, :);
    ms{col} = m_;
end

end