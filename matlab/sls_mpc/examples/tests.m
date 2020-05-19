clear all; close all; clc;

%% Algorithm 1 setup plant + parameters
rng(420);
sys1    = LTISystem;
sys1.Nx = 8; 
alpha = 0.8; rho = 1; actDens = 0.6; 
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

params1.Q_ = eye(sys1.Nx);
params1.R_ = eye(sys1.Nu);

%% Algorithm 1, no constraints
params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Alg1, no constraints');

%% Algorithm 1, with state ub
params1.stateConsMtx_ = eye(sys1.Nx);
params1.stateUB_      = 1;

params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Alg1, state ub');

%% Algorithm 1, with state lb + ub
params1.stateConsMtx_ = eye(sys1.Nx);
params1.stateUB_      = 1;
params1.stateLB_      = 0;

params1.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys1, x01, params1);

params1.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys1, x01, params1);

printAndPlot(params1, x, u, xVal, uVal, 'Alg1, state lb + ub');

%% Algorithm 2 setup plant + parameters
sys2    = LTISystem;
sys2.Nx = 4; 
alpha = 0.8; rho = 1; actDens = 0.5; 
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
params2.Q_ = diag(ones(Nx,1)) + diag(-1/2*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
params2.R_ = eye(sys2.Nu);

% Constraints
for i = 1:2:2*(Nx-1)
    K(i,i)     = 1; 
    K(i,i+2)   = -1;
    K(i+1,i)   = -1; 
    K(i+1,i+2) = 1;
end
K            = K(1:Nx,1:Nx); 
K(Nx-1:Nx,:) = zeros(2,Nx);

%% Algorithm 2, no constraints
params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Alg2, no constriants');

%% Algorithm 2, with state ub
params2.stateConsMtx_ = K;
params2.stateUB_      = 0.8;

params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Alg2, with state ub');

%% Algorithm 2, with state lb + ub
params2.stateConsMtx_ = K;
params2.stateUB_      = 0.8;
params2.stateLB_      = -0.8;

params2.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x02, params2);

params2.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x02, params2);

printAndPlot(params2, x, u, xVal, uVal, 'Alg2, with state lb + ub');
 
%% Local function to print values + plot graphs
function printAndPlot(params, x, u, xVal, uVal, myTitle)
    % Calculate costs + plot 
    tSim   = params.tHorizon_;
    obj    = get_cost_fn(params, x, u);
    objVal = get_cost_fn(params, xVal, uVal);

    % Print costs (sanity check: should be close)
    fprintf('Distributed cost: %f\n', obj);
    fprintf('Centralized cost: %f\n', objVal);

    figure()
    plot(1:tSim+1,xVal(1,:),'b',1:tSim+1,x(1,:),'*b')
    xlabel('Time');
    ylabel('State (1st only)');
    legend('Centralized', 'Distributed');
    title(myTitle);
end