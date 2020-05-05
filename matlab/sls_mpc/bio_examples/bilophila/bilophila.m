clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
q = 0.8;
k = [1; 1; 1];

% initial conditions
x0 = [1e-3; 1e-3; 1e-3; 1; 1e-3];

% simulation length
tHorizon = 20;

% SLS length
tFIR = 4;

% disturbance (step)
val1 = 1e-3;
dur1 = 8;
val2 = 2;
ws   = [val1 * ones(1, dur1), val2 * ones(1, tHorizon-dur1)];

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup system
Nx = 5;
Nu = 3;

sys = LTISystem();
sys.Nx = Nx; sys.Nu = Nu;

xs      = zeros(Nx, tHorizon);
us      = zeros(Nu, tHorizon);
xs(:,1) = x0;

u_lb = -1 * ones(Nu, 1);
u_ub = [0; 1; 0];
x_lb = 1e-3 * ones(Nx, 1);

for t=1:tHorizon-1
    fprintf('Time: %d\n', t);
    
    x_ = xs(:,t);
    u_ = us(:,t);
    w_ = ws(t);

    [Ac, Bc]        = linearize_bilo(x_, u_, q, k);    
    [sys.A, sys.B2] = discretize(Ac, Bc, Ts);
    
    % coordinate shift bounds
    u_tilde_lb = u_lb - u_;
    u_tilde_ub = u_ub - u_;
    x_tilde_lb = x_lb - x_;
    
    y_shift = pinv(sys.A) * f_bilo(x_, u_, w_, q, k)*Ts;
    yt      = y_shift;
    y_lb    = x_tilde_lb + y_shift;
    
    [y, u_tilde] = mpc_bilo(sys, tFIR, u_tilde_lb, u_tilde_ub, y_lb, yt);
    
    % calculate x(t+1) from dynamics and u_tilde directly
    x_tilde_nxt = f_bilo(x_, u_, w_, q, k)*Ts + sys.B2*u_tilde;
    xs(:,t+1) = x_ + x_tilde_nxt;
    
    % update actuation
    us(:,t) = u_ + u_tilde;    
end

%% Plot
figure(); 
hold on;

for i=1:5
    plot(1:tHorizon, xs(i,:));
end

legend('x1 (taurine)', 'x2 (sulfite)', 'x3 (sulfide)', 'x4 (energy)', 'x5 (growth)'); 
