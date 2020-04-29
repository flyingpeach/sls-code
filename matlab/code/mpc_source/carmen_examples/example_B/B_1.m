% Algorithm I, System B
clear all; close all; clc;

%% Setup plant + sweep variables
setup_plant_b;

rng(2020);
x0 = rand(sys.Nx, 1);

localities = [3 5 7 10];
numLocs    = length(localities);
times      = zeros(1, numLocs);
timeCents  = zeros(1, numLocs);
iters      = zeros(1, numLocs);

%% MPC parameters
params = MPCParams();

params.tFIR_     = 5;
params.tHorizon_ = 10;
params.maxIters_ = 5000;
params.rho_      = 5;
params.eps_p_    = 1e-4;
params.eps_d_    = 1e-3;
params.solnMode_ = MPCSolMode.ClosedForm;

params.Q_ = eye(sys.Nx);
params.R_ = eye(sys.Nu);

%% Sweep over locality sizes
for i=1:numLocs
    params.locality_ = localities(i);
    
    % Distributed MPC
    [x, u, times(i), iters(i)] = mpc_algorithm_1(sys, x0, params);
    
    % Centralized MPC (for validation + comparison)
    [xVal, uVal, timeCents(i)] = mpc_centralized(sys, x0, params);   
end

%% Plot
figure(1)
subplot(1,2,1)
plot(localities, times,'m-s','LineWidth',2)
hold on
plot(localities, timeCents,'b-s','LineWidth',2)
xlabel('$$\#\ pendulums\ in\ the\ network$$','Interpreter','latex','Fontsize', 10)
ylabel('$$Avg\ runtime\ per\ MPC\ iteration\ for\ each\ state\ (s)$$','Interpreter','latex','Fontsize', 10)
leg4 = legend('$$Localized\ ADMM\ Solution$$', '$$Centralized\ Solution$$');
set(leg4,'Interpreter','latex','Fontsize', 8);

subplot(1,2,2)
plot(localities, iters,'m-s','LineWidth',2)
xlabel('$$\#\ pendulums\ in\ the\ network$$','Interpreter','latex','Fontsize', 10)
ylabel('$$Avg\ \#\ ADMM\ iters\ per\ MPC\ iteration\ for\ each\ state\ (s)$$','Interpreter','latex','Fontsize', 10)
