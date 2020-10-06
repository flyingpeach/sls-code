function [clMaps1, clMaps2] = two_blend_sls(sys, slsParams, satParams)
% Solve for two linear SLS closed-loop maps as per
% https://arxiv.org/pdf/2006.12766.pdf    

T   = slsParams.T_;
eta = satParams.eta_;

cvx_begin quiet

if isempty(slsParams.constraints_)
    variable R1s(sys.Nx, sys.Nx, T)
    variable R2s(sys.Nx, sys.Nx, T)
    variable M1s(sys.Nu, sys.Nx, T)
    variable M2s(sys.Nu, sys.Nx, T)
else
    expression R1s(sys.Nx, sys.Nx, T)
    expression R2s(sys.Nx, sys.Nx, T)
    expression M1s(sys.Nu, sys.Nx, T)
    expression M2s(sys.Nu, sys.Nx, T)
end
    
% populate decision variables for ease-of-use
R1 = cell(T, 1); R2 = cell(T, 1);
M1 = cell(T, 1); M2 = cell(T, 1);
for t=1:T
    R1{t} = R1s(:,:,t); R2{t} = R2s(:,:,t);
    M1{t} = M1s(:,:,t); M2{t} = M2s(:,:,t);
end

if ~isempty(slsParams.constraints_)
    [RSupp, MSupp, count] = get_supports(sys, slsParams);
    variable RMSuppVals1(count)
    variable RMSuppVals2(count)

    [R1, M1] = add_sparse_constraints(R1, M1, RSupp, MSupp, RMSuppVals1, T);
    [R2, M2] = add_sparse_constraints(R2, M2, RSupp, MSupp, RMSuppVals2, T);
end

% achievability constraints
R1{1} == eye(sys.Nx); R2{1} == eye(sys.Nx);
R1{T} == zeros(sys.Nx, sys.Nx);
R2{T} == zeros(sys.Nx, sys.Nx);
for t=1:T-1
    R1{t+1} == sys.A*R1{t} + sys.B2*M1{t};
    R2{t+1} == sys.A*R2{t} + sys.B2*M2{t};
end

R1Concat = []; R2Concat = [];
M1Concat = []; M2Concat = [];
B1Blk = [];
for t=1:T
    R1Concat = [R1Concat R1{t}];  R2Concat = [R2Concat R2{t}];
    M1Concat = [M1Concat M1{t}];  M2Concat = [M2Concat M2{t}];
    B1Blk = blkdiag(B1Blk, sys.B1);
end

% get covariance matrix from noise (sigma_w in paper)
eyeSize = sys.Nx * slsParams.T_;
covMtx  = get_cov_mtx(satParams, eyeSize);

% objective
RM = [R1Concat  R2Concat;
      M1Concat  M2Concat];
objectiveMtx = [sys.C1, sys.D12] * RM * sqrtm(covMtx);
objective    = norm(objectiveMtx, 'fro');

% saturation constraints
R1Inf = norm(R1Concat*B1Blk, Inf);  R2Inf = norm(R2Concat*B1Blk, Inf);
M1Inf = norm(M1Concat*B1Blk, Inf);  M2Inf = norm(M2Concat*B1Blk, Inf);

if ~isinf(satParams.xMax_)
    eta(1)*R1Inf + (eta(2) - eta(1))*R2Inf <= satParams.xMax_;
end
if ~isinf(satParams.uMax_)
    eta(1)*M1Inf + (eta(2) - eta(1))*M2Inf <= satParams.uMax_;
end

minimize(objective)
cvx_end

clMaps1              = CLMaps();
clMaps1.R_           = R1;
clMaps1.M_           = M1;
clMaps1.solveStatus_ = cvx_status;

clMaps2              = CLMaps();
clMaps2.R_           = R2;
clMaps2.M_           = M2;
clMaps2.solveStatus_ = cvx_status; 

if strcmp(cvx_status, 'Solved')
    fprintf(['Solved!', '\n\n']);
else
    sls_warning(['Solver exited with status ', cvx_status]);
end
