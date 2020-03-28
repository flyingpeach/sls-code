function [F, G] = get_ctrller_constraint(sys, clMaps, Tc)
% Returns the matrices F, G that provide the constraint F[Rc{2:Tc}; Mc] = G
% (This assumes Rc{1} = I, as expected)
% Inputs
%    sys    : LTISystem containing system matrices
%    clMaps : contains closed-loop maps (R, M)
%    Tc     : order of controller matrices

T      = length(clMaps.R_);
Rblock = cell(T+Tc, Tc+1);
Mblock = cell(T+Tc, Tc+1);

% initialize with zeros
for i=1:T+Tc 
    for j=1:Tc+1
        Rblock{i,j} = zeros(sys.Nx, sys.Nx);
        Mblock{i,j} = zeros(sys.Nu, sys.Nx);
    end
end

for j=1:Tc+1
    for idx=1:T
        Rblock{idx+j-1,j} = clMaps.R_{idx};
        Mblock{idx+j-1,j} = clMaps.M_{idx};
    end
end

Ralpha = zeros(Tc*sys.Nx, (Tc+1)*sys.Nx);
Rbeta  = zeros(T *sys.Nx, (Tc+1)*sys.Nx);
Malpha = zeros(Tc*sys.Nu, (Tc+1)*sys.Nx);
Mbeta  = zeros(T *sys.Nu, (Tc+1)*sys.Nx);

for j=1:Tc+1
    jx = get_range(j, sys.Nx);    
    for i=1:Tc %"alpha" blocks are first Tc rows
        ix = get_range(i, sys.Nx);
        iu = get_range(i, sys.Nu);        
        Ralpha(ix, jx) = Rblock{i,j};
        Malpha(iu, jx) = Mblock{i,j};
    end
    for i=1:T %"beta" blocks are last T rows
        ix = get_range(i, sys.Nx);
        iu = get_range(i, sys.Nu);        
        Rbeta(ix, jx) = Rblock{Tc+i,j};
        Mbeta(iu, jx) = Mblock{Tc+i,j};
    end
end

Dellc = zeros((Tc+1)*sys.Nx, Tc*(sys.Nx + sys.Nu));
for idx=1:Tc
    idx1 = get_range(idx, sys.Nx);
    idx2 = get_range(idx+1, sys.Nx);
    idx3 = get_range(idx+1, sys.Nx);
    idx4 = Tc * sys.Nx + get_range(idx, sys.Nu);
    
    Dellc(idx1, idx1) = eye(sys.Nx);
    Dellc(idx2, idx1) = -sys.A;
    
    Dellc(idx3, idx4) = -sys.B2;
end

myEye = eye(Tc*(sys.Nx + sys.Nu));
myZer = zeros(T*(sys.Nx + sys.Nu), Tc*(sys.Nx + sys.Nu));

FG = [Ralpha; Malpha; Rbeta; Mbeta] * Dellc - [myEye; myZer];
F  = FG(:,sys.Nx+1:end);
G  = -FG(:, 1:sys.Nx);