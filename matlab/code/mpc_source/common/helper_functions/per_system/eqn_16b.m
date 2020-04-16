function Psi_loc = eqn_16b(phi_coli, lamb_coli, ZABi, Eyei)
           
M       = ZABi'*pinv(ZABi*ZABi');
Psi_loc = (phi_coli + lamb_coli) + M*(Eyei - ZABi*(phi_coli + lamb_coli));

end