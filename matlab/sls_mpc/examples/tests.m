clear all; close all; clc;

%% Algorithm 1 setup plant + parameters
rng(420);
sys1    = LTISystem;
sys1.Nx = 8; 
alpha = 0.8; rho = 2; actDens = 0.6; 
generate_dbl_stoch_chain(sys1, rho, actDens, alpha);
x01 = rand(sys1.Nx, 1);

params1 = MPCParams();
params1.locality_ = 3;
params1.tFIR_     = 5;
params1.tHorizon_ = 10;
params1.maxIters_ = 5000;
params1.rho_      = 1; 
params1.eps_p_    = 1e-3;
params1.eps_d_    = 1e-3;

params1.QSqrt_ = diag(randi([1, 3], sys1.Nx, 1));
params1.RSqrt_ = diag(randi([1, 3], sys1.Nu, 1));

plotStates = 1;
plotInputs = 2;

%% TEST A: Algorithm 1, no constraints
params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Test A', plotStates, plotInputs);

%% TEST B: Algorithm 1, with state constraints
params1.stateConsMtx_ = eye(sys1.Nx);
params1.stateUB_      = 1.7;  % not a tight constraint
params1.stateLB_      = -0.5; % tight constraint

params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Test B', plotStates, plotInputs);

%% TEST C: Algorithm 1, with state+input constraints
% state constraints still apply from TEST B if run sequentially
params1.inputConsMtx_ = eye(sys1.Nu);
params1.inputLB_      = -1.5; % tight constraint

params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Test C', plotStates, plotInputs);

%% Algorithm 2 setup plant + parameters
rng(420);
sys2    = LTISystem;
sys2.Nx = 4; 
alpha = 0.8; rho = 2; actDens = 0.7; 
generate_dbl_stoch_chain(sys2, rho, actDens, alpha);

x02 = rand(sys2.Nx, 1);

params2 = MPCParams();
params2.locality_ = 3;
params2.tFIR_     = 5;
params2.tHorizon_ = 10;
params2.maxIters_ = 5000;
params2.rho_      = 1; 
params2.eps_p_    = 1e-3;
params2.eps_d_    = 1e-3;

params2.maxItersCons_ = 500;
params2.mu_           = 1;
params2.eps_x_        = 1e-3;
params2.eps_z_        = 1e-3;

Nx = sys2.Nx;
params2.QSqrt_ = diag(ones(Nx,1)) + diag(-1/2*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
params2.RSqrt_ = diag(randi([1, 3], sys2.Nu, 1));

% State constraints (different coupling from cost)
K = zeros(Nx);
K(1,:) = [1 1 0 0];
K(2,:) = [0 1 1 0];

plotStates = [2 3];
plotInputs = 1:sys2.Nu;

%% TEST D: Algorithm 2, no constraints
params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Test D', plotStates, plotInputs);

%% TEST E: Algorithm 2, with state constraints
params2.stateConsMtx_ = K;
params2.stateUB_      = 0.8; % tight constraint

params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Test E', plotStates, plotInputs);

%% TEST F: Algorithm 2, with state+input constraints
% state constraints still apply from TEST E if run sequentially
params2.inputConsMtx_ = eye(sys2.Nu);
params2.inputLB_      = -0.8; % tight constraint

params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Test F', plotStates, plotInputs);
 
%% Local function to print values + plot graphs
function printAndPlot(params, x, u, xVal, uVal, myTitle, states, inputs)
    % Calculate costs + plot 
    tSim   = params.tHorizon_;
    obj    = get_cost_fn(params, x, u);
    objVal = get_cost_fn(params, xVal, uVal);

    % Print costs (sanity check: should be close)
    fprintf('Distributed cost: %f\n', obj);
    fprintf('Centralized cost: %f\n', objVal);
    
    time = 1:tSim;
    figure();
    subplot(2,1,1); hold on;

    title(myTitle);
    for i=1:length(states)
        state = states(i);
        plot(time, xVal(state,time),'b', time, x(state,time),'*b');
    end
    stateLabel = ['states ', num2str(states)];
    ylabel(stateLabel);
    legend('Centralized', 'Distributed')
    
    subplot(2,1,2); hold on;
    for i=1:length(inputs)
        input = inputs(i);
        plot(time, uVal(input,:),'g', time, u(input,:),'*g');
    end
    inputLabel = ['inputs ', num2str(inputs)];
    ylabel(inputLabel);
    xlabel('Time');
end