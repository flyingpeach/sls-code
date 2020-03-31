function lqrCost = get_ctrller_lqr_cost(sys, ctrller, tTotal)
% Result will be incorrect if ctrller not stable or tTotal too small

[R_, M_]     = get_cl_map(sys, ctrller, tTotal);
objectiveMtx = get_objective_mtx(sys, R_, M_);
lqrCost      = get_H2_norm(objectiveMtx);
