clear all; close all; clc;

%% Setup plant + parameters
rng(420);
Nx     = 8; alpha = 0.8; rho = 2; actDens = 0.6; 
sys    = generate_dbl_stoch_chain(Nx, rho, actDens, alpha);
sys.B1 = eye(sys.Nx);

tHorizon = 10;
x0       = zeros(sys.Nx, 1);

paramsNom = MPCParams();
paramsNom.locality_ = 3;
paramsNom.tFIR_     = 10;
paramsNom.maxIters_ = 5000;
paramsNom.rho_      = 1; 
paramsNom.eps_p_    = 1e-2;
paramsNom.eps_d_    = 1e-2;

paramsNom.QSqrt_ = diag(ones(Nx,1)) + diag(-1/2*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
paramsNom.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));

% Disturbance
dist = [-1 1 1 -1 -1 -1 1 1 1 1];
w    = ones(sys.Nx, 1) * dist;

%% Case A: Unconstrained nominal (as sanity check)
paramsNom.mode_        = MPCMode.Centralized;
[xsCentA, usCentA, ~] = sls_mpc(sys, x0, w, paramsNom, tHorizon);

%% Case B(1): Constrained centralized robust
paramsRob = copy(paramsNom);

% Bounded state
KState = zeros(sys.Nx);
KState(1, 1) = 1; % constraint on x1 + x2
KState(1, 2) = 1;

paramsRob.stateConsMtx_ = KState;
paramsRob.stateUB_ = 3 * ones(sys.Nx, 1); % tight constraint
paramsRob.stateLB_ = -3.5 * ones(sys.Nx, 1); % tight constraint

% Bounded disturbance
paramsRob.distConsMtx_ = eye(sys.Nx);
paramsRob.distUB_      =  1 * ones(sys.Nx, 1);
paramsRob.distLB_      = -1 * ones(sys.Nx, 1);

% Robust Centralized
paramsRob.mode_       = MPCMode.Centralized;
[xsCentB, usCentB, ~] = sls_mpc(sys, x0, w, paramsRob, tHorizon);

%% Case B(2): Constrained distributed robust
% Adaptive ADMM
paramsRob.tau_i_   = 2;
paramsRob.tau_d_   = 2;
paramsRob.muAdapt_ = 10;
paramsRob.rhoMax_  = 80;

paramsRob.mode_    = MPCMode.Distributed;
[xsB, usB, statsB] = sls_mpc(sys, x0, w, paramsRob, tHorizon);

fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsB.time_, statsB.iters_);

%% Postprocessing
objCentA = get_cost_fn(paramsRob, xsCentA, usCentA);
objCentB = get_cost_fn(paramsRob, xsCentB, usCentB);
objB     = get_cost_fn(paramsRob, xsB, usB);

fprintf('Unconstrained cost (cent): %f\n',    objCentA);
fprintf('Constrained cost   (cent): %f\n',    objCentB);
fprintf('Constrained cost   (dist): %f\n',    objB);

plotState = 1;
plotInput = 2;

time = 1:size(xsCentA, 2);
figure();
subplot(2,1,1); hold on;
plot(time, xsCentA(1, time) + xsCentA(2, time), 'k');
plot(time, xsCentB(1, time) + xsCentB(2, time), 'b');
plot(time, xsB(1, time) + xsB(2, time), '*b');
stateLabel = ['State 1 + 2'];
ylabel(stateLabel);

subplot(2,1,2); hold on;
plot(time, usCentA(plotInput, time), 'k');
plot(time, usCentB(plotInput, time), 'g');
plot(time, usB(plotInput, time), '*g');

inputLabel = ['Input ', num2str(plotInput)];
ylabel(inputLabel);
xlabel('Time');
legend('Unconstr (Cent)', 'Constr (Cent)', 'Constr (Dist)');
