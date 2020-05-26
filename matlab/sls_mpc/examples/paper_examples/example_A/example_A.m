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
params.eps_p_    = 1e-4;
params.eps_d_    = 1e-3;
params.RSqrt_    = eye(Nu);

% Used in algorithm 2 only
params.eps_x_    = 1e-3;
params.eps_z_    = 1e-4;

tHorizon = 60;

%% Open Loop
xOpen(:,1) = x0;
for t = 1:tHorizon-1
    xOpen(:,t+1) = sys.A*xOpen(:,t);
end
figure(1)
plot(1:tHorizon, xOpen(1,:), 1:tHorizon, xOpen(3,:));

%% Case 1
params.maxIters_ = 5000;
params.rho_      = 5;
params.QSqrt_    = eye(Nx);

params.mode_        = MPCMode.Distributed;
[x1, u1, ~]         = sls_mpc(sys, x0, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCent1, uCent1, ~] = sls_mpc(sys, x0, params, tHorizon);

plot_a(params, x1, u1, xCent1, uCent1, 'Case 1');

%% Case 2
params.maxIters_ = 5000;
params.rho_      = 1;

params.maxItersCons_ = 200;
params.mu_           = 1;

params.QSqrt_ = zeros(Nx,Nx);
for i = 1:2:Nx
    params.QSqrt_(i,i) = 1;
    if i > 1
        params.QSqrt_(i,i-2) = -1/2;
    end
    if i < Nx-1
        params.QSqrt_(i,i+2) = -1/2;
    end
    if i < Nx
        params.QSqrt_(i+1,i+1) = .01;
    end
end

params.mode_        = MPCMode.Distributed;
[x2, u2, ~]         = sls_mpc(sys, x0, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCent2, uCent2, ~] = sls_mpc(sys, x0, params, tHorizon);

plot_a(params, x2, u2, xCent2, uCent2, 'Case 2');

%% Case 3
params.maxIters_ = 5000;
params.rho_      = 300;

params.maxItersCons_ = 500;
params.mu_           = 50;

for i = 1:2:2*(Nx-1)
    K(i,i)     = 1; 
    K(i,i+2)   = -1;
    K(i+1,i)   = -1; 
    K(i+1,i+2) = 1;
end
K            = K(1:Nx,1:Nx); 
K(Nx-1:Nx,:) = zeros(2,Nx);

params.stateConsMtx_ = K;
params.stateUB_      = 0.05;

params.mode_        = MPCMode.Distributed;
[x3, u3, ~]         = sls_mpc(sys, x0, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCent3, uCent3, ~] = sls_mpc(sys, x0, params, tHorizon);

plot_a(params, x3, u3, xCent3, uCent3, 'Case 3');
