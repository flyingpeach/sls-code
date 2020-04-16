% Algorithm II, System B
clear all; close all; clc;

localities = [3 5 7 10];
numLocs    = length(localities);
times      = zeros(1, numLocs);
timeCents  = zeros(1, numLocs);
iters      = zeros(1, numLocs);

numPendula = 10;

% MPC parameters
rho   = 1; 
eps_d = 1e-3; % convergence criterion for ||Psi(k+1) - Psi(k)||
eps_p = 1e-4; % convergence criterion for ||Phi(k+1) - Psi(k+1)||
maxIters = 10000;

mu       = 1;
eps_pres = 1e-3;
eps_dres = 1e-4;
maxConsensusIters = 500;

for index=1:numLocs
    locality = localities(index);
    clear x u x_VAL u_VAL LocalityR LocalityM
    
    setup_system_b;
    
    % Weights
    Q = diag(ones(Nx,1))+diag(-1/2*ones(Nx-2,1),2)+diag(-1/2*ones(Nx-2,1),-2);
    S = eye(Nu);

    % Cost matrix
    C = [];
    for t = 0:tFIR-1
        C = blkdiag(C, Q);
    end
    for t = 0:tFIR-2
        C = blkdiag(C, S);
    end    

    % Coupling
    indices = cell(1, length(C));
    for i = 1:length(C)
        for j = 1:length(C)
            if C(i,j) ~= 0
                indices{i} = [indices{i} j];
            end
        end
    end
    
    % Distributed MPC   
    [x, u, times(index), iters(index)] = mpc_algorithm2(Nx, Nu, A, B, d, tFIR, tSim, x0, ...
                              eps_d, eps_p, rho, maxIters, ...
                              indices, eps_pres, eps_dres, mu, ...
                              maxConsensusIters, C);

    % Centralized MPC (for validation + comparison)
    [xVal, uVal, timeCents(idx)] = mpc_centralized(Nx, Nu, A, B, d, ...
        Q, S, tFIR, tSim, x0);
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
