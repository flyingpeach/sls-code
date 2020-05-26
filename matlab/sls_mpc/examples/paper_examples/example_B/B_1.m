% Algorithm I, System B
clear all; close all; clc;
rng(2020);

%% MPC parameters
params = MPCParams();

params.tFIR_     = 5;
params.maxIters_ = 5000;
params.rho_      = 5;
params.eps_p_    = 1e-4;
params.eps_d_    = 1e-3;

tHorizon = 10;

%% Sweep over locality sizes
numPendula     = 10;
sys            = setup_plant_b(numPendula);
x0             = rand(sys.Nx, 1);
params.QSqrt_  = eye(sys.Nx);
params.RSqrt_  = eye(sys.Nu);

localities   = [3 5 7 10];
numLocs      = length(localities);
times_l      = zeros(1, numLocs);
timeCents_l  = zeros(1, numLocs);
iters_l      = zeros(1, numLocs);

for i=1:numLocs
    params.locality_ = localities(i);
    
    params.mode_ = MPCMode.Distributed;
    [x, u, times_l(i), iters_l(i)] = sls_mpc(sys, x0, params, tHorizon);

    params.mode_ = MPCMode.Centralized;
    [xVal, uVal, timeCents_l(i)] = sls_mpc(sys, x0, params, tHorizon);   
end

%% Plot sweep over locality sizes
localityBool = 1;
plot_b(localities, times_l, timeCents_l, iters_l, localityBool);

%% Sweep over network sizes
params.locality_ = 3;

sizes       = [10 50 100 200];
numSizes    = length(sizes);
times_n     = zeros(1, numSizes);
timeCents_n = zeros(1, numSizes);
iters_n     = zeros(1, numSizes);

% index 1 of both sweeps have exactly the same parameters; don't rerun mpc
times_n(1)     = times_l(1);
iters_n(1)     = iters_l(1);
timeCents_n(1) = timeCents_l(1);

for i=2:numSizes
    numPendula = sizes(i);
    sys        = setup_plant_b(numPendula);
    x0         = rand(sys.Nx, 1);
    
    params.QSqrt_  = eye(sys.Nx);
    params.RSqrt_  = eye(sys.Nu);

    params.mode_ = MPCMode.Distributed;
    [x, u, times_n(i), iters_n(i)] = sls_mpc(sys, x0, params, tHorizon);

    params.mode_ = MPCMode.Centralized;
    [xVal, uVal, timeCents_n(i)] = sls_mpc(sys, x0, params, tHorizon);
end

%% Plot sweep over network sizes
localityBool = 0;
plot_b(sizes, times_n, timeCents_n, iters_n, localityBool);
