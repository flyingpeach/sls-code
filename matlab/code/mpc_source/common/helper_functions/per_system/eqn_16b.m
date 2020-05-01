function psi_ = eqn_16b(phi_, lamb_, zab_, eye_)
% phi_  : column of phi
% lamb_ : column of Lambda
% zab_  : column of ZAB (part of sls constraints)
% eye_  : column of shifted identity matrix (part of sls constraints)
           
M    = zab_'*pinv(zab_*zab_');
psi_ = (phi_ + lamb_) + M*(eye_ - zab_*(phi_ + lamb_));

end