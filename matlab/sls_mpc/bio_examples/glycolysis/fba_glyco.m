%% Flux Balance Analysis
clear all; clc
%% Parameters

Nx = 3;

% model parameters
q = 0.8;
k = [2 4 0.1];
c = [0 0 -1]';

% initial conditions
xs(:,1) = [0.6; 0.6; 0.1];

% simulation length
tHorizon = 10;

% disturbance (step)
val1 = 0;  % value of first disturbance
dur1 = 5;  % how long first disturbance lasts
val2 = -3; % value of second disturbance (lasts forever)
ws = [val1*ones(1, dur1) val2*ones(1, tHorizon-dur1)];

% stoichiometry matrix
S = [1 -1 0; -q q+1 -1; 0 0 1];

% sampling time for discretization
Ts = 0.1;

%% LP

for t = 1:tHorizon-1

    cvx_begin

    variable v(Nx,1)

    minimize (c'*v)
    S*[k(1)*v(1);k(2)*v(2);k(3)*v(3)] == -[0 1 0]'*ws(t);
    zeros(Nx,1) <= v;
    v <= exp(ones(Nx,1));

    cvx_end
    
    flux = [k(1)*v(1);k(2)*v(2);k(3)*v(3)];
    
    %% Compute the state

    deltax = (S*flux+[0 1 0]'*ws(t+1))*Ts;
    xs(:,t+1) = xs(:,t)+deltax;
    
end
%%
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
