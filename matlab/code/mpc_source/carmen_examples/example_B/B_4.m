% Algorithm II, System B
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
params.rho_      = 300;
params.eps_p_    = 1e-4;
params.eps_d_    = 1e-3;
params.solnMode_ = MPCSolMode.UseSolver;

params.maxItersCons_ = 500;
params.mu_           = 50;
params.eps_x_        = 1e-3;
params.eps_z_        = 1e-4;

params.constrUpperbnd_ = 0.5;

params.Q_ = diag(ones(Nx,1)) + diag(-1/2*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
params.R_ = eye(sys.Nu);

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

%% Sweep over locality sizes
for i=1:numLocs
    params.locality_ = localities(i);

    % Distributed MPC
    [x, u, times(i), iters(i)] = mpc_algorithm_2(sys, x0, params);

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
