function spectralRadius = check_int_stab(sys, ctrller)

intStabMtx     = get_int_stab_mtx(sys, ctrller);
spectralRadius = max(abs(eig(intStabMtx))); % biggest eigenvalue
