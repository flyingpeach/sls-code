clear all; close all; clc;

%% Setup plant + parameters
rng(400);
sys    = LTISystem;
sys.Nx = 10; alpha = 0.8; rho = 1.5; actDens = 0.5; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);
sys.B1 = eye(sys.Nx);

tHorizon = 10;
x0       = rand(sys.Nx, 1);
dist = [-1 1 1 -1 -1 -1 1 1 1 1];
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

%% Synthesize w/o terminal set
fprintf('Doing rMPC with *no* terminal set...\n')
xCentA = nan(sys.Nx, tHorizon);
uCentA = nan(sys.Nu, tHorizon);
xCentA(:,1) = x0;

for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t+1, tHorizon);
    [~, uCentA(:,t), ~] = mpc_centralized(sys, xCentA(:,t), params);
    
    if any(isnan(uCentA(:,t)))
        mpc_warning('rMPC solver failed/infeasible!');
        break;
    end
    
    xCentA(:,t+1) = sys.A*xCentA(:,t) + sys.B2*uCentA(:,t) + sys.B1*w(:,t);
end

%% Add terminal set + synthesize w/ terminal set
fprintf('Synthesizing robust terminal set...\n')
[params, tStats] = terminal_set(sys, params);
fprintf('Robust terminal set: avgTime: %.4f, avgIters: %.4f\n\n', tStats.time_, tStats.iters_);

fprintf('Doing rMPC *with* terminal set...\n')
params.mode_        = MPCMode.Centralized;
[xCentB, uCentB, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Distributed;
[xB, uB, statsB]    = sls_mpc(sys, x0, w, params, tHorizon);

fprintf('rMPC w/terminal set: avgTime: %.4f, avgIters: %.4f\n\n', statsB.time_, statsB.iters_);

objTermDist = get_cost_fn(params, xB, uB);
objTermCent = get_cost_fn(params, xCentB, uCentB);
fprintf('Distributed cost w/terminal set: %f\n', objTermDist);
fprintf('Centralized cost w/terminal set: %f\n', objTermCent);

%% Plotting
plotState = 2;
plotInput = 1;

time = 1:tHorizon;
figure();
subplot(2,1,1); hold on;
plot(time, xCentB(plotState, :), 'b');
plot(time, xB(plotState, :), '*b');
plot(time, xCentA(plotState, :), 'r');
ylabel('State');
legend('Centralized [Term]', 'Distributed [Term]', 'Centralized [NoTerm]')
    
subplot(2,1,2); hold on;
plot(time, uCentB(plotInput, :), 'g');
plot(time, uB(plotInput, :), '*g');
plot(time, uCentA(plotInput, :), 'm');
ylabel('Input');
legend('Centralized [Term]', 'Distributed [Term]', 'Centralized [NoTerm]')
xlabel('Time');