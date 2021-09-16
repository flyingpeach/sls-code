function [G, g] = get_constraint_g(params)
% This function works with time-invariant constraints of the form
%    distLB <= distConsMtx*disturbance <= distUB
% and converts them into the form of G * delta <= g

T = params.tFIR_;

G_half = [];
g_ub   = []; % portion of g representing upper bounds
g_lb   = []; % portion of g representing lower bounds
for t = 1:T-1
    G_half = blkdiag(G_half, params.distConsMtx_);
    g_ub = [g_ub; params.distUB_];
    g_lb = [g_lb; params.distLB_];
end

G = [G_half; -G_half];
g = [g_ub; -g_lb];

end
