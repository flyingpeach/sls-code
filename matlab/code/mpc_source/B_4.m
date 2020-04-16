% Algorithm II, System B
clear all; close all; clc;

localities = [3 5 7 10];
numLocs    = length(localities);
times      = zeros(1, numLocs);
timeCents  = zeros(1, numLocs);
iters      = zeros(1, numLocs);

numPendula = 10;
% MPC parameters
rho   = 300;
eps_d = 1e-3; % convergence criterion for ||Psi(k+1) - Psi(k)||
eps_p = 1e-4; % convergence criterion for ||Phi(k+1) - Psi(k+1)||
maxIters = 5000;

mu       = 50;
eps_pres = 1e-3;
eps_dres = 1e-4;
maxConsensusIters = 500;
up = .5;

for index=1:numLocs
    locality = localities(index);
    clear x u x_VAL u_VAL LocalityR LocalityM
    
    setup_system_b;

    % Weights
    Q = diag(ones(Nx,1))+diag(-1/2*ones(Nx-2,1),2)+diag(-1/2*ones(Nx-2,1),-2);
    S = diag(ones(Nu,1));
        
    % Cost matrix
    C = [];
    for t = 0:tFIR-1
        C = blkdiag(C, Q);
    end
    for t = 0:tFIR-2
        C = blkdiag(C, S);
    end    

    % Constraints
    for i = 1:2:2*(Nx-1)
        K1(i,i) = 1; K1(i,i+2) = -1;
        K1(i+1,i) = -1; K1(i+1,i+2) = 1;
    end
    K1 = K1(1:Nx,1:Nx); K1(Nx-1:Nx,:) = zeros(2,Nx);

    Ksmall = zeros(2*Nx,Nx); j = 0;
    for i = 1:2*Nx
        if mod(i,4) == 1 | mod(i,4) == 2
            j = j + 1;
            Ksmall(i,:) = K1(j,:);
        else
        end
    end

    % Build the big constraint matrix
    K = [zeros(size(Ksmall)) zeros(2*Nx,(tFIR-1)*Nx)];
    K =  [K; zeros(2*Nx,Nx) zeros(size(Ksmall)) zeros(2*Nx,(tFIR-2)*Nx)];
    for t = 2:tFIR-1
        K =  [K; zeros(2*Nx,t*Nx) Ksmall zeros(2*Nx,(tFIR-t-1)*Nx)];
    end

    % Coupling (since the coupling from cost is more extensive than the coupling from the constraints this
                % time we only count the coupling from the cost)
    indices = cell(1, length(C));
    for i = 1:length(C)
        for j = 1:length(C)
            if C(i,j) ~= 0
                indices{i} = [indices{i} j];
            end
        end
    end
    
    % Distributed MPC   
    [x, u, times(index), iters(index)] = mpc_algorithm2_gurobi(Nx, Nu, A, B, d, tFIR, tSim, x0, ... % system
                              eps_d, eps_p, rho, maxIters, ...
                              indices, eps_pres, eps_dres, mu, ...
                              maxConsensusIters, C, up, K);
    
    % Centralized MPC (for validation + comparison)
    coupling = true;
    [xVal, uVal, timesCents(idx)] = mpc_centralized(Nx, Nu, A, B, d, Q, S, tFIR, tSim, x0, up, Ksmall, coupling);
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
