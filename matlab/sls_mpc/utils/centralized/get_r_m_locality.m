function [r_loc, m_loc] = get_r_m_locality(sys, locality)
% r_loc and m_loc contain sparsity patterns on R, M dictated by locality
% in the locality-only case, r_loc{k} is the same for all k

Comms_Adj = abs(sys.A)>0;
r_loc     = Comms_Adj^(locality-1) > 0;
m_loc     = abs(sys.B2)'*r_loc > 0;

end