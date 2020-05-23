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

%% Algorithm 1, no constraints
params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Alg1, no constraints', 1);

%% Algorithm 1, with state lb
params1.stateConsMtx_ = eye(sys1.Nx);
params1.stateLB_      = -0.5;

params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Alg1, state lb', 1);

%% Algorithm 1, with state lb + ub
params1.stateConsMtx_ = eye(sys1.Nx);
params1.stateUB_      = 1.7;
params1.stateLB_      = -0.5;

params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Alg1, state lb + ub', 1);

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

% Constraints
K = zeros(Nx);
K(1,1) = 1; K(1,3) = 1;
K(2,2) = 1; K(2,4) = 1;

%% Algorithm 2, no constraints
params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Alg2, no constraints', [1,3]);

%% Algorithm 2, with state ub
params2.stateConsMtx_ = K;
params2.stateUB_      = 0.6;

params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Alg2, with state ub', [1,3]);

%% Algorithm 2, with state lb + ub
params2.stateConsMtx_ = K;
params2.stateUB_      = 0.6;
params2.stateLB_      = 0;

params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Alg2, with state lb + ub', [1,3]);
 
%% Local function to print values + plot graphs
function printAndPlot(params, x, u, xVal, uVal, myTitle, states)
    % Calculate costs + plot 
    tSim   = params.tHorizon_;
    obj    = get_cost_fn(params, x, u);
    objVal = get_cost_fn(params, xVal, uVal);

    % Print costs (sanity check: should be close)
    fprintf('Distributed cost: %f\n', obj);
    fprintf('Centralized cost: %f\n', objVal);
    
    figure(); hold on;
    
    for i=1:length(states)
        state = states(i);
        plot(1:tSim+1,xVal(state,:),'b',1:tSim+1,x(state,:),'*b')
    end
    
    xlabel('Time');
    ylabel('States');
    legend('Centralized', 'Distributed');
    title(myTitle);
end