function IZA_ZB = get_sls_constraints(sys, tFIR)

Nx = sys.Nx; Nu = sys.Nu;

I = kron(eye(tFIR),eye(Nx));

Z = kron(eye(tFIR-1),eye(Nx));
Z = [zeros(Nx,Nx*(tFIR));Z,zeros(Nx*(tFIR-1),Nx)];

A_repmat = repmat({sys.A},tFIR,1);
B_repmat = repmat({sys.B2},tFIR,1);

A_block  = blkdiag(A_repmat{:});
B_block  = blkdiag(B_repmat{:});

IZA_ZB = [I - Z*A_block, -Z*B_block];
IZA_ZB = IZA_ZB(:,1:end-Nu);
end