clear all; close all; clc;

%% Setup plant + parameters
Nx     = 10; alpha = 0.8; rho = 1.5; actDens = 0.6; 
sys    = generate_dbl_stoch_chain(Nx, rho, actDens, alpha);
sys.B1 = eye(sys.Nx);

tHorizon = 12;
x0   = 0.2*ones(sys.Nx, 1);
dist = [-1 1 1 -1 -1 -1 1 1 1 1 -1 -1 -1 -1 -1];
w    = 0.1 * ones(sys.Nx, 1) * dist;

% Params
params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 2;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 1e-3;
params.eps_d_    = 1e-3;

params.distConsMtx_ = eye(sys.Nx);
params.distUB_      =  1 * ones(sys.Nx, 1);
params.distLB_      = -1 * ones(sys.Nx, 1);

params.stateConsMtx_     = eye(sys.Nx);
params.stateUB_          = 5 * ones(sys.Nx, 1);
params.stateLB_          = -params.stateUB_;    
params.QSqrt_ = eye(sys.Nx);
params.RSqrt_ = eye(sys.Nu);

% Consensus (for terminal cost only)
params.mu_           = 1.2;
params.eps_x_        = 1e-3;
params.eps_z_        = 1e-3;
params.maxItersCons_ = 1000;

%% Case A: No terminal set
fprintf('Doing MPC with NO terminal set...\n')
xCentA = nan(sys.Nx, tHorizon);
uCentA = nan(sys.Nu, tHorizon);
xCentA(:,1) = x0;

for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t+1, tHorizon);
    [~, uCentA(:,t), ~] = mpc_centralized(sys, xCentA(:,t), params);
    
    if any(isnan(uCentA(:,t)))
        mpc_warning('MPC solver failed/infeasible!');
        break;
    end
    
    xCentA(:,t+1) = sys.A*xCentA(:,t) + sys.B2*uCentA(:,t) + sys.B1*w(:,t);
end

%% Case B: Terminal constraint only
fprintf('Synthesizing terminal set...\n')
[params, tStats] = terminal_set(sys, params);
fprintf('Terminal set: avgTime: %.4f, avgIters: %.4f\n\n', tStats.time_, tStats.iters_);

fprintf('Doing MPC with terminal constraint...\n')
params.mode_        = MPCMode.Centralized;
[xCentB, uCentB, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Distributed;
[xB, uB, statsB]    = sls_mpc(sys, x0, w, params, tHorizon);

fprintf('MPC, terminal constraint\n');
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsB.time_, statsB.iters_);

objDistB = get_cost_fn(params, xB, uB);
objCentB = get_cost_fn(params, xCentB, uCentB);
fprintf('Dist. cost w/terminal constraint: %f\n', objDistB);
fprintf('Cent. cost w/terminal constraint: %f\n', objCentB);

%% Case C: Terminal constraint + cost
params.terminal_cost_ = true;

fprintf('Doing MPC with terminal constraint + cost...\n')
params.mode_        = MPCMode.Centralized;
[xCentC, uCentC, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Distributed;
[xC, uC, statsC]    = sls_mpc(sys, x0, w, params, tHorizon);
 
fprintf('MPC, terminal constraint+cost:\n');
fprintf('avgTime: %.4f, avgIters: %.4f, avgConsIters: %.4f\n\n', statsC.time_, statsC.iters_, statsC.consIters_);

objDistC = get_cost_fn(params, xC, uC);
objCentC = get_cost_fn(params, xCentC, uCentC);
fprintf('Dist. cost w/terminal constraint + cost: %f\n', objDistC);
fprintf('Cent. cost w/terminal constraint + cost: %f\n', objCentC);

%% Plotting
plotState = 1;
plotInput = 1;

time = 1:tHorizon;
figure();
subplot(2,1,1); hold on;
plot(time, xCentA(plotState, :), 'r');
plot(time, xCentB(plotState, :), 'b');
plot(time, xB(plotState, :), '*b');
plot(time, xCentC(plotState, :), 'g');
plot(time, xC(plotState, :), '*g');
ylabel('State');
legend('Cent [CaseA]', 'Cent [CaseB]', 'Dist [CaseB]', 'Cent [CaseC]', 'Dist [CaseC]');

subplot(2,1,2); hold on;
plot(time, uCentA(plotInput, :), 'r');
plot(time, uCentB(plotInput, :), 'b');
plot(time, uB(plotInput, :), '*b');
plot(time, uCentC(plotInput, :), 'g');
plot(time, uC(plotInput, :), '*g');
ylabel('Input');
xlabel('Time');
