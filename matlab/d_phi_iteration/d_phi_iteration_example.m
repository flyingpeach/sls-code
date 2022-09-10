clear; clc; close;

% Plant parameters
Nx      = 10;
Nu      = Nx;
rho     = 3.0; % Spectral radius

% Generate random bidirectional ring
A = zeros(Nx, Nx); rng(420);
A = A + diag(rand(Nx-1,1)-0.5,1);
A = A + diag(rand(Nx-1,1)-0.5,-1);
%A(1,Nx) = rand()-0.5; A(Nx, 1) = rand()-0.5; % bidirectional ring
A(1,Nx) = 1; A(Nx, 1) = 1; % bidirectional ring
A = A / max(abs(eig(A))) * rho; % give A the desired spectral radius

sys.A  = A; 
sys.B2 = eye(Nx);

params = DPhiParams();
params.locality_ = 3; % communicate with neighbors and neighbors of neighbors
params.tFIR_     = 30;
params.betaStep_ = 5e-2; % minimum robust bound progress per step

% Regulated output z = x + B*u, corresponds to per-node self uncertainty
params.Hx_       = eye(Nx);
params.Hu_       = sys.B2;

% LQR penalty
params.QSqrt_ = eye(Nx);
params.RSqrt_ = sqrt(50)*eye(Nu);

%% Get LQR values to compare against
A = sys.A; B = sys.B2;
Q = params.QSqrt_'*params.QSqrt_;
R = params.RSqrt_'*params.RSqrt_;
S = idare(A, B, Q, R);
K = -(B'*S*B + R)\(B'*S*A); % by convention MATLAB uses (A-BK)

Nx = size(A, 1);

costLQR = 0;
for i = 1:Nx
    ei = zeros(Nx,1); ei(i) = 1;
    costLQR = costLQR + ei'*S*ei;
end

Phi_xLQR{1} = eye(Nx); Phi_uLQR{1} = K;
MAX_HORIZON = 1000; % Hardcoded, assumes horizon < MAX_HORIZON
for k=1:MAX_HORIZON
    Phi_xLQR{k+1} = A*Phi_xLQR{k} + B*Phi_uLQR{k};
    Phi_uLQR{k+1} = K*Phi_xLQR{k+1};
end

M  = make_m(params, Phi_xLQR, Phi_uLQR);
[~, betaLQRNu]   = d_step_minimize(M, StabType.Nu);
[~, betaLQRL1]   = d_step_minimize(M, StabType.L1);
[~, betaLQRLInf] = d_step_minimize(M, StabType.LInf);

%% Nu Scan
params.stabType_     = StabType.Nu;

% Alg1
params.randomizeD_ = false;
[Phi_xsNu1, Phi_usNu1, betasNu1, costsNu1] = d_phi_iteration(sys, params);

% Alg2
params.randomizeD_ = true;
[Phi_xsNu2, Phi_usNu2, betasNu2, costsNu2] = d_phi_iteration(sys, params);

%% L1 Scan
params.stabType_   = StabType.L1;

% Alg1
params.randomizeD_ = false;
[Phi_xsL1_1, Phi_usL1_1, betasL1_1, costsL1_1] = d_phi_iteration(sys, params);

% Alg2
params.randomizeD_ = true;
[Phi_xsL1_2, Phi_usL1_2, betasL1_2, costsL1_2] = d_phi_iteration(sys, params);
 
%% LInf Scan
params.stabType_   = StabType.LInf;

% Alg1
params.randomizeD_ = false;
[Phi_xsLInf_1, Phi_usLInf_1, betasLInf_1, costsLInf_1] = d_phi_iteration(sys, params);

% Alg2
params.randomizeD_ = true;
[Phi_xsLInf_2, Phi_usLInf_2, betasLInf_2, costsLInf_2] = d_phi_iteration(sys, params);

%% Plot everything
figure(1); 

ylims = [1.17 1.31];

idxNu1 = 1:18; idxNu2 = 1:18;
subplot(1,3,1); hold on;
plot(betaLQRNu./betasNu1(idxNu1), costsNu1(idxNu1)/costLQR, '*-');
plot(betaLQRNu./betasNu2(idxNu2), costsNu2(idxNu2)/costLQR, 'o-');
legend('Nu, Alg1', 'Nu, Alg2'); 
ylabel('Normalized Cost');
xlabel('Normalized Robustness Margin')
ylim(ylims)

idxL1_1 = 1:30; idxL1_2 = 5:35;
subplot(1,3,2); hold on;
plot(betaLQRL1./betasL1_1(idxL1_1), costsL1_1(idxL1_1)/costLQR, '*-');
plot(betaLQRL1./betasL1_2(idxL1_2), costsL1_2(idxL1_2)/costLQR, 'o-');
legend('L1, Alg1', 'L1, Alg2'); 
xlabel('Normalized Robustness Margin')
ylim(ylims)

idxLInf_1 = 1:25; idxLInf_2 = 8:35;
subplot(1,3,3); hold on;
plot(betaLQRLInf./betasLInf_1(idxLInf_1), costsLInf_1(idxLInf_1)/costLQR, '*-');
plot(betaLQRLInf./betasLInf_2(idxLInf_2), costsLInf_2(idxLInf_2)/costLQR, 'o-');
legend('LInf, Alg1', 'LInf, Alg2'); 
xlabel('Normalized Robustness Margin')
ylim(ylims)