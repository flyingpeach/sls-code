% Algorithm I, System A
clear all; close all; clc;

%% Setup plant
setup_plant_a;

rng(2020);
x0 = rand(sys.Nx, 1);

%% MPC Parameters
params = MPCParams();

params.locality_ = 3;
params.tFIR_     = 20;
params.tHorizon_ = 30;
params.maxIters_ = 5000;
params.rho_      = 5;
params.eps_p_    = 1e-4;
params.eps_d_    = 1e-3;
params.solnMode_ = MPCSolMode.ClosedForm;

params.Q_ = eye(Nx);
params.R_ = eye(Nu);

%% Distributed MPC
[x, u, ~] = mpc_algorithm_1(sys, x0, params);

%% Centralized MPC (for validation + comparison)
[xVal, uVal, ~] = mpc_centralized(sys, x0, params);

%% Calculate costs + plot 
tSim   = params.tHorizon_;
obj    = get_cost_fn(params, x, u);
objVal = get_cost_fn(params, xVal, uVal);

% Output costs
fprintf('Distributed cost: %f\n', obj);
fprintf('Centralized cost: %f\n', objVal);

figure(1)
plot(1:tSim+1,xVal(1,:),'b',1:tSim+1,x(1,:),'*b',1:tSim+1,xVal(3,:),'g',1:tSim+1,x(3,:),'*g')
xlabel('$$Time$$','interpreter','latex','Fontsize', 10)
ylabel('$$\theta_{1},\ \theta_{2}$$','Interpreter','Latex','Fontsize', 10)
leg1 = legend('$$\theta_{1}\ Centralized\ MPC$$', '$$\theta_{1}\ Localized\ MPC\ using\ ADMM$$','$$\theta_{2}\ Centralized\ MPC$$', '$$\theta_{2}\ Localized\ MPC\ using\ ADMM$$');
set(leg1,'Interpreter','latex'); set(leg1, 'Fontsize', 8)
title('Subsystems 1 and 2')