clear all; close all; clc;

%% Setup plant
rng(420);

sys    = LTISystem;
sys.Nx = 4; 
alpha  = 0.8; rho = 1; actDens = 1; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

x0 = rand(sys.Nx, 1);

%% MPC Parameters
params = MPCParams();

params.locality_ = 3;
params.tFIR_     = 20;
params.tHorizon_ = 10;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 1e-3;
params.eps_d_    = 1e-3;
params.solnMode_ = MPCSolMode.UseSolver;

params.maxItersCons_ = 500;
params.mu_           = 1;
params.eps_x_        = 1e-3;
params.eps_z_        = 1e-3;

Nx = sys.Nx;
params.Q_ = diag(ones(Nx,1)) + diag(-1/2*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
params.R_ = eye(sys.Nu);

params.constrUpperbnd_ = 0.5;

% Constraints
for i = 1:2:2*(Nx-1)
    K1(i,i)     = 1; 
    K1(i,i+2)   = -1;
    K1(i+1,i)   = -1; 
    K1(i+1,i+2) = 1;
end
K1            = K1(1:Nx,1:Nx); 
K1(Nx-1:Nx,:) = zeros(2,Nx);
  
Ksmall = zeros(2*Nx,Nx); j = 0;
for i = 1:2*Nx
    if mod(i,4) == 1 || mod(i,4) == 2
        j = j + 1;
        Ksmall(i,:) = K1(j,:);
    else
    end
end

params.constrMtx_ = Ksmall;

[x1, u1, ~]       = mpc_algorithm_2(sys, x0, params);

%%
[xVal1, uVal1, ~] = mpc_centralized(sys, x0, params);

printAndPlot(params, x1, u1, xVal1, uVal1, 'Alg1 with Solver');

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
