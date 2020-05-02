function [lqrCost, intRad, l1norm] = get_ctrller_stats(sys, ctrller, tTotal)
% Result will be incorrect if ctrller not stable or tTotal too small

lqrCost = get_ctrller_lqr_cost(sys, ctrller, tTotal);
intRad  = get_int_stab_radius(sys, ctrller);
objMtx  = get_objective_mtx(sys, ctrller.Rc_, ctrller.Mc_);
l1norm  = get_L1_norm(objMtx);
end