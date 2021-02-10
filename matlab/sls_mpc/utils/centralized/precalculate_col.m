function [zabs, eyes, zabis] = precalculate_col(sys, T, s_c)

Nx = sys.Nx;

Eye = [eye(Nx); zeros(Nx*(T-1),Nx)];
ZAB = get_sls_constraint(sys, T);

zabs  = cell(Nx, 1);
eyes  = cell(Nx, 1);
zabis = cell(Nx, 1);

for col = 1:Nx % ith column
    zab_       = ZAB(:, s_c{col});
    zeroRows   = find(all(zab_ == 0, 2));
    keepRows   = setdiff(1:size(ZAB, 1), zeroRows);
    zabs{col}  = ZAB(keepRows, s_c{col});
    eyes{col}  = Eye(keepRows, col);
    zabis{col} = zabs{col}'*pinv(zabs{col}*zabs{col}');
end

end