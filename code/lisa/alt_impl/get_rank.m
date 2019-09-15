function [rankF, rankF2] = get_rank(sys, slsParams, slsOuts, Tc, tol)
% Returns the rank of F and F2. This can be used to determine whether an
% alternate implementation exists, and how big the solution space is
% Outputs
%    rankF     : rank of F
%    rankF2    : rank of F without its first column
% Inputs
%    sys       : LTISystem containing system matrices
%    slsParams : SLSParams containing parameters
%    slsOuts   : contains info from SLS (original R, M)
%    Tc        : length of the approximate solution
%    tol       : tolerance to use for rank calculation

T      = slsParams.tFIR_;
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
        Rblock{idx+j-1,j} = slsOuts.R_{idx};
        Mblock{idx+j-1,j} = slsOuts.M_{idx};
    end
end

Ralpha = zeros(Tc*sys.Nx, (Tc+1)*sys.Nx);
Rbeta  = zeros(T *sys.Nx, (Tc+1)*sys.Nx);
Malpha = zeros(Tc*sys.Nu, (Tc+1)*sys.Nx);
Mbeta  = zeros(T *sys.Nu, (Tc+1)*sys.Nx);

% helper function to convert cells of block matrices into giant matrix
get_range = @(idx, size) (size*(idx-1)+1:size*(idx-1)+size);

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

F = [Ralpha; Malpha; Rbeta; Mbeta] * Dellc - [myEye; myZer];
rankF  = rank(F, tol);
rankF2 = rank(F(:,sys.Nx+1:end), tol);
