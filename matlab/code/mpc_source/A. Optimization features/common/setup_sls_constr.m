E1 = [eye(Nx);zeros(Nx*(tFIR-1),Nx)];

I = kron(eye(tFIR),eye(Nx));

Z = kron(eye(tFIR-1),eye(Nx));
Z = [zeros(Nx,Nx*(tFIR));Z,zeros(Nx*(tFIR-1),Nx)];

tmp = repmat({A},tFIR,1);
Aa = blkdiag(tmp{:});
clear tmp
tmp = repmat({B},tFIR,1);
Bb = blkdiag(tmp{:});

IZAa = I - Z*Aa;
ZB = -Z*Bb;

IZA_ZB = [IZAa ZB];
IZA_ZB = IZA_ZB(:,1:end-Nu);