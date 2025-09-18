function pPhix = pI_m_A(Phix,SLSPoles,A)
n = size(Phix{1},1);
Polenum = size(Phix,1);
cvx_begin
expression pPhix(n,n*Polenum)
for l = 1:Polenum
    pPhix(:,(l-1)*n+1:l*n) = (SLSPoles(l)*eye(n)-A)*Phix{l};
end
cvx_end
end