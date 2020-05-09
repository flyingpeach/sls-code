clear all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
q = 0.8;
k = [2 4 0.1];

% initial conditions
x0 = [0.6; 0.6; 0.1];

% simulation length
tHorizon = 5;

% disturbance (step)
val1 = 0;  % value of first disturbance
dur1 = 2;  % how long first disturbance lasts
val2 = -3; % value of second disturbance (lasts forever)
ws = [val1*ones(1, dur1) val2*ones(1, tHorizon-dur1)];

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
% stoichiometry matrix
Nx = 3;
S = [ 1   -1  0;
     -q  q+1 -1];

K = diag(k);

xs      = zeros(Nx, tHorizon);
fluxes  = zeros(Nx, tHorizon);
xs(:,1) = x0;

e2 = [0 1]';

for t = 1:tHorizon-1

    cvx_begin

    variable v(Nx,1)
    maximize([0 0 1] * v);
    subject to
    S*K*v + e2*ws(t) == 0;

    v >= zeros(Nx, 1);
    v <= exp(ones(Nx, 1));

    cvx_end
    
    flux = K*v;
    
    % update state
    deltax      = (S*flux + e2*ws(t+1)) * Ts;
    xs(1:2,t+1) = xs(1:2,t) + deltax;
    xs(3,t+1)   = xs(3,t) + flux(3)*Ts;
end
%% Plot
figure();

subplot(2,1,1);
hold on;
for i=1:2
    plot(1:tHorizon, xs(i,:));
end
ylabel('x1, x2');
legend('x1 (intermeds)', 'x2 (atp)');

subplot(2,1,2);
plot(1:tHorizon, xs(3,:));
ylabel('x3 (growth)');
