% Need MPT3 for this script
clear; clc;

Nx = 8;
Nu = 3;

A = 4*rand(Nx, Nx);
B = rand(Nx, Nu);

Q = eye(Nx);
R = 1e-2*eye(Nu);
K = dlqr(A, B, Q, R);

A_CL = A - B*K;
model = LTISystem('A', A_CL);

% Constraints on states
model.x.min = -5 * ones(Nx, 1);
model.x.max = 5 * ones(Nx, 1);

% Compute max positive invariant set using mpt3 (this is correct)
trueSet = model.invariantSet();

% Use our method
params = MPCParams();
params.QSqrt_ = Q;
params.RSqrt_ = sqrt(R);
params.locality_ = 10; % global

params.stateConsMtx_ = eye(Nx);
params.stateUB_      = 5 * ones(Nx, 1);
params.stateLB_      = -params.stateUB_;

[H_term, h_term] = terminal_set(A, B, Nx, params);

ourSet = Polyhedron('A', H_term, 'b', h_term);
ourSet.minHRep(); 

% Check that our sets are equal
norm(sort(trueSet.H) - sort(ourSet.H))


% Look at sets
% trueSet.H
% [H_term h_term]








