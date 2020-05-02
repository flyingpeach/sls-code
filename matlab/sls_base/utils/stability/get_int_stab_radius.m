function spectralRadius = get_int_stab_radius(sys, ctrller)

intStabMtx     = get_int_stab_mtx(sys, ctrller);
spectralRadius = max(abs(eig(intStabMtx))); % biggest eigenvalue
