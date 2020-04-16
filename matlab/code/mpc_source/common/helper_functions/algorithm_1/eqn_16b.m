function Psi_loc = eqn_16b(ci, s_ci, phi_coli, lamb_coli, ZAB, Eye)

ZABi     = ZAB(:, sci);
zeroRow  = find(all(ZABi == 0, 2));
keepIdxs = setdiff(linspace(1,Nx*tFIR,Nx*tFIR), zeroRow);
ZABi     = ZAB(keepIdxs, s_ci); 
Eyei     = Eye(keepIdxs, ci);
            
M       = ZABi'*pinv(ZABi*ZABi');
Psi_loc = (phi_coli + lamb_coli) + M*(Eyei - ZABi*(phi_coli + lamb_coli));

end