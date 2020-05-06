clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
q = 0.8;
k = 1;
a = 1;

% initial conditions
x0 = [1.5; 1.5];
z0 = log(x0);

% simulation length
tHorizon = 20;

% SLS length
tFIR = 4;

% disturbance (step)
val1 = 1e-2;
dur1 = 10;
val2 = 1e-2;
ws   = [val1 * ones(1, dur1), val2 * ones(1, tHorizon-dur1)];

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup system
Nz = 2;
Nu = 2;

sys = LTISystem();
sys.Nx = Nz; sys.Nu = Nu;

zs      = zeros(Nz, tHorizon);
us      = zeros(Nu, tHorizon);
zs(:,1) = z0;

u_ub = zeros(Nu, 1);
x_lb = 1e-3 * ones(Nz, 1);
z_lb = log(x_lb);

for t=1:tHorizon-1
    fprintf('Time: %d\n', t);
    
    z_ = zs(:,t);
    u_ = us(:,t);
    w_ = ws(t);

    [Ac, Bc]        = linearize_glyco_log(z_, u_, q, k, a);    
    [sys.A, sys.B2] = discretize(Ac, Bc, Ts);
    
    % coordinate shift bounds
    u_tilde_ub = u_ub - u_;
    z_tilde_lb = z_lb - z_;
    
    y_shift = pinv(sys.A - eye(Nz)) * f_glyco_log(z_, u_, w_, q, k, a)*Ts;
    yt      = y_shift;
    y_lb    = z_tilde_lb + y_shift;
    
    [y, u_tilde] = mpc_glyco(sys, tFIR, y_lb, u_tilde_ub, yt);
    
    % calculate x(t+1) from dynamics and u_tilde directly
    z_tilde_nxt = f_glyco_log(z_, u_, w_, q, k, a)*Ts + sys.B2*u_tilde;
    zs(:,t+1)   = z_ + z_tilde_nxt;
    
    % update actuation
    us(:,t) = u_ + u_tilde;    
end

%% Plot
figure(); 
hold on;

xs = exp(zs);

for i=1:2
    plot(1:tHorizon, xs(i,:));
end

legend('x1 (intermeds)', 'x2 (atp)');
