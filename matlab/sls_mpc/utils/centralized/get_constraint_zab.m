function ZAB = get_constraint_zab(sys, T)

Nx = sys.Nx; Nu = sys.Nu;

I = kron(eye(T),eye(Nx));

Z = kron(eye(T-1),eye(Nx));
Z = [zeros(Nx,Nx*(T));Z,zeros(Nx*(T-1),Nx)];

A_repmat = repmat({sys.A},T,1);
B_repmat = repmat({sys.B2},T,1);

A_block  = blkdiag(A_repmat{:});
B_block  = blkdiag(B_repmat{:});

ZAB = [I - Z*A_block, -Z*B_block];
ZAB = ZAB(:,1:end-Nu);
end