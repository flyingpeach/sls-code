clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameter
q = 0.8;

% initial conditions
x0 = 0.1*ones(5, 1);

% simulation length
tHorizon = 10;

% SLS length
tFIR = 5;

% disturbance
ws = [0.1 * ones(1, tHorizon/2), ones(1, tHorizon/2)];

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup system
u_lb  = [-1; -1; -1];
u_ub  = [0; 1; 0];
x4_lb = 1e-3;
    
Nx = 5;
Nu = 3;

sys = LTISystem();
sys.Nx = Nx; sys.Nu = Nu;

xs      = zeros(Nx, tHorizon);
us      = zeros(Nu, tHorizon);
xs(:,1) = x0;

e1 = [1 0 0 0 0]';

for t=1:tHorizon-1
    fprintf('Time: %d\n', t);
    
    x_ = xs(:,t);
    u_ = us(:,t); % always zero
    w_ = ws(:,t);

    [Ac, Bc]        = linearize_bilo(x_, u_, q);    
    [sys.A, sys.B2] = discretize(Ac, Bc, Ts);
    
    % coordinate shift bounds
    u_tilde_lb  = u_lb - u_;
    u_tilde_ub  = u_ub - u_;
    x4_tilde_lb = x4_lb - x_(4);
    
    % y, u_tilde are coordinate shifts of x, u
    y0 = pinv(sys.A) * f_bilo(x_, u_, w_, q);
    
    [y, u_tilde] = mpc_bilo(sys, tFIR, u_tilde_lb, u_tilde_ub, x4_tilde_lb, y0);

    % calculate x(t+1) from dynamics and u_tilde directly
    x_tilde_nxt = f_bilo(x_, u_, w_, q)*Ts + sys.B2*u_tilde + e1*(ws(:,t+1) - w_);
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
