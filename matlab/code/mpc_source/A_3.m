% Algorithm II, System A
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
params.rho_      = 300;
params.eps_p_    = 1e-4;
params.eps_d_    = 1e-3;
params.solnMode_ = MPCSolMode.UseSolver;

params.maxItersCons_ = 500;
params.mu_           = 50;
params.eps_x_        = 1e-3;
params.eps_z_        = 1e-4;

params.stateUpperbnd_ = 0.05;
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
params.R_ = eye(Nu);

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

params.constraintMtx_ = Ksmall;

%% Distributed MPC
[x, u, ~] = mpc_algorithm2(sys, x0, params);

%% Validation
coupling = true;
[xVal, uVal, ~] = mpc_centralized(Nx, Nu, A, B, d, Q_, S, tFIR, tSim, x0, up, KSmall, coupling);

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