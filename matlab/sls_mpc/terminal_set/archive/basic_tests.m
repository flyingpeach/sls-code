% Many called functions use the mpt3 toolbox! Please install it
clear all; close all; clc;

%% Check fourier_motzkin
% Example from wiki
A = [2 -5 4;
     3 -6 3;
     -1 5 -2;
     -3 2 6];

b = [10 9 -7 12]';

[A_tilde, b_tilde] = fourier_motzkin(A, b);

%% One-step controllable set
% Borelli textbook example 10.2 (p.188)
sys    = LTISystem();
sys.A  = [1.5 0;
          1   -1.5];
sys.B2 = [1; 0];

sys.Nx = 2; sys.Nu = 1;

params = MPCParams();

params.stateConsMtx_ = eye(sys.Nx);
params.stateUB_      =  10 * ones(sys.Nx, 1);
params.stateLB_      = -10 * ones(sys.Nx, 1); 

params.inputConsMtx_ = 1;
params.inputUB_      = 5;
params.inputLB_      = -5;

N = 1;
[H1_, h1_] = nominal_ctrb_set(sys, params, N);

%% Max control invariant set
% Borelli textbook example 10.6 (p.195)
[Hmax_, hmax_, iters] = nominal_max_ci_set(sys, params);

[Hmax_, hmax_]
iters



