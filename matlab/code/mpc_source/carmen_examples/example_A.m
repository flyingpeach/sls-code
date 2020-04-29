clear all; close all; clc;

%% Setup plant
sys = setup_plant_a();
Nu = sys.Nu; Nx = sys.Nx; % for reduced notation

rng(2020);
x0 = rand(sys.Nx, 1);

%% MPC Parameters (common)
params = MPCParams();

params.locality_ = 3;
params.tFIR_     = 20;
params.tHorizon_ = 60;
params.eps_p_    = 1e-4;
params.eps_d_    = 1e-3;

% Used in algorithm 2 only
params.eps_x_        = 1e-3;
params.eps_z_        = 1e-4;

params.R_ = eye(Nu);

% For ease of notation
tSim = params.tHorizon_;

%% Open Loop
xOpen(:,1) = x0;
for t = 1:tSim
    xOpen(:,t+1) = sys.A*xOpen(:,t);
end
figure(1)
plot(1:tSim+1, xOpen(1,:), 1:tSim+1, xOpen(3,:));

%% Case 1
params.maxIters_ = 5000;
params.rho_      = 5;
params.solnMode_ = MPCSolMode.ClosedForm;

params.Q_ = eye(Nx);

[x1, u1, ~]       = mpc_algorithm_1(sys, x0, params);
[xVal1, uVal1, ~] = mpc_centralized(sys, x0, params);

printAndPlot(params, x1, u1, xVal1, uVal1, 'Case 1');

%% Case 2
params.maxIters_ = 5000;
params.rho_      = 1;
params.solnMode_ = MPCSolMode.ClosedForm;

params.maxItersCons_ = 200;
params.mu_           = 1;

params.Q_ = zeros(Nx,Nx);
for i = 1:2:Nx
    params.Q_(i,i) = 1;
    if i > 1
        params.Q_(i,i-2) = -1/2;
    end
    if i < Nx-1
        params.Q_(i,i+2) = -1/2;
    end
    if i < Nx
        params.Q_(i+1,i+1) = .01;
    end
end

[x2, u2, ~]       = mpc_algorithm_2(sys, x0, params);
[xVal2, uVal2, ~] = mpc_centralized(sys, x0, params);

printAndPlot(params, x2, u2, xVal2, uVal2, 'Case 2');

%% Case 3
params.maxIters_ = 5000;
params.rho_      = 300;
params.solnMode_ = MPCSolMode.UseSolver;

params.maxItersCons_ = 500;
params.mu_           = 50;

params.stateUpperbnd_ = 0.05;

% Constraints
for i = 1:2:2*(Nx-1)
    K1(i,i)     = 1; 
    K1(i,i+2)   = -1;
    K1(i+1,i)   = -1; 
    K1(i+1,i+2) = 1;
end
K1            = K1(1:Nx,1:Nx); 
K1(Nx-1:Nx,:) = zeros(2,Nx);
  
params.constraintMtx_ = zeros(2*Nx,Nx); j = 0;
for i = 1:2*Nx
    if mod(i,4) == 1 || mod(i,4) == 2
        j = j + 1;
        params.constraintMtx_(i,:) = K1(j,:);
    else
    end
end

[x3, u3, ~]       = mpc_algorithm_2(sys, x0, params);
[xVal3, uVal3, ~] = mpc_centralized(sys, x0, params);

printAndPlot(params, x3, u3, xVal3, uVal3, 'Case 3');

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
    plot(1:tSim+1,xVal(1,:),'b',1:tSim+1,x(1,:),'*b',1:tSim+1,xVal(3,:),'g',1:tSim+1,x(3,:),'*g')
    xlabel('$$Time$$','interpreter','latex','Fontsize', 10)
    ylabel('$$\theta_{1},\ \theta_{2}$$','Interpreter','Latex','Fontsize', 10)
    leg1 = legend('$$\theta_{1}\ Centralized\ MPC$$', '$$\theta_{1}\ Localized\ MPC\ using\ ADMM$$','$$\theta_{2}\ Centralized\ MPC$$', '$$\theta_{2}\ Localized\ MPC\ using\ ADMM$$');
    set(leg1,'Interpreter','latex'); set(leg1, 'Fontsize', 8)
    title(strcat(myTitle, ', Subsystems 1 and 2'));
end
