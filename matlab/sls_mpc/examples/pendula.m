clear all; close all; clc;

%% Setup plant
n  = 4; % number of pendulums
Nx = 2*n; 
Nu = n;

m = 1; k = 1; f = 3; g = 10; l = 1;
block_off_diag  = [0    0; k*l/m  f/(m*l)];
block_diag      = [0 1; -g-2*k*l/m -2*f/(m*l)];

Ac = zeros(Nx,Nx); j = 0; % A matrix (continuous time)
for i = 1:2:Nx
    j = j+1;
    if j == 1 % first node
        Ac (i:i+1,i+2:i+3) = block_off_diag;
        Ac (i:i+1,i:i+1) = block_diag;
    elseif j == Nx/2  % last node      
        Ac (i:i+1,i:i+1) = block_diag;
        Ac (i:i+1,i-2:i-1) = block_off_diag;
    else
        Ac (i:i+1,i+2:i+3) = block_off_diag;
        Ac (i:i+1,i:i+1) = block_diag;
        Ac (i:i+1,i-2:i-1) = block_off_diag;
    end
end

% B matrix (continous time)
Bc = zeros(Nx, Nu); j = 0;
for i = 1:2:Nx
    j = j+1;
    Bc (i:i+1,j) = [0; 1];
end

% Discretize + set up system
Ts = .1;

sys     = LTISystem();
sys.Nx  = Nx;
sys.Nu  = Nu;
sys.A   = (eye(Nx)+Ac*Ts);
sys.B2  = Ts*Bc;

sys.sanity_check_mpc();

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

tHorizon = 50;

%% Open Loop
xOpen(:,1) = x0;
for t = 1:tHorizon-1
    xOpen(:,t+1) = sys.A*xOpen(:,t);
end
figure(1)
plot(1:tHorizon, xOpen(1,:), 1:tHorizon, xOpen(3,:));

%% Case 1 (no coupling)
params.maxIters_ = 5000;
params.rho_      = 5;
params.QSqrt_    = eye(Nx);

params.mode_        = MPCMode.Distributed;
[x1, u1, ~]         = sls_mpc(sys, x0, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCent1, uCent1, ~] = sls_mpc(sys, x0, params, tHorizon);

plot_pendula(params, x1, u1, xCent1, uCent1, 'Case 1');

%% Case 2 (coupled objective)
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

plot_pendula(params, x2, u2, xCent2, uCent2, 'Case 2');

%% Local function to assist with plotting
function plot_pendula(params, x, u, xVal, uVal, myTitle)

% Calculate costs + plot 
time   = 1:size(x,2);
obj    = get_cost_fn(params, x, u);
objVal = get_cost_fn(params, xVal, uVal);

% Print costs (sanity check: should be close)
fprintf('Distributed cost: %f\n', obj);
fprintf('Centralized cost: %f\n', objVal);

figure()
plot(time,xVal(1,:),'b',time,x(1,:),'*b',time,xVal(3,:),'g',time,x(3,:),'*g')
xlabel('$$Time$$','interpreter','latex','Fontsize', 10)
ylabel('$$\theta_{1},\ \theta_{2}$$','Interpreter','Latex','Fontsize', 10)
leg1 = legend('$$\theta_{1}\ Centralized\ MPC$$', '$$\theta_{1}\ Localized\ MPC\ using\ ADMM$$','$$\theta_{2}\ Centralized\ MPC$$', '$$\theta_{2}\ Localized\ MPC\ using\ ADMM$$');
set(leg1,'Interpreter','latex'); set(leg1, 'Fontsize', 8)
title(strcat(myTitle, ', Subsystems 1 and 2'));

end
