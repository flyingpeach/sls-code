function psi_ = eqn_16b(phi_, lamb_, zab_, eye_, zabi_)
% phi_  : column of phi
% lamb_ : column of Lambda
% zab_  : column of ZAB (part of sls constraints)
% eye_  : column of shifted identity matrix (part of sls constraints)
% zabi_ : pseudoinverse of zab

psi_ = (phi_ + lamb_) + zabi_*(eye_ - zab_*(phi_ + lamb_));

end