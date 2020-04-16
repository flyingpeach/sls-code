% Algorithm I, System B
clear all; close all; clc;

localities = [3 5 7 10];
numLocs    = length(localities);
times      = zeros(1, numLocs);
timeCents  = zeros(1, numLocs);
iters      = zeros(1, numLocs);

rho   = 10; % admm parameter
eps_d = 1e-3; % convergence criterion for ||Psi(k+1) - Psi(k)||
eps_p = 1e-4; % convergence criterion for ||Phi(k+1) - Psi(k+1)||
maxIters = 10000;

% bounds for Gurobi solver
up  =  1.2;
low = -0.2;

numPendula = 3;

for idx=1:numLocs
    locality = localities(idx);
    setup_system_b;
    Q = eye(Nx);
    S = eye(Nu);
    % Distributed MPC
    [x, u, times(idx), iters(idx)] = mpc_algorithm1_gurobi(Nx, Nu, A, B, d, tFIR, tSim, x0, ...
                 eps_d, eps_p, rho, maxIters, up, low);
    
    % Centralized MPC (for validation + comparison)
    [xVal, uVal, timeCents(idx)] = mpc_centralized(Nx, Nu, A, B, d, ...
        Q, S, tFIR, tSim, x0, up, low);
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
