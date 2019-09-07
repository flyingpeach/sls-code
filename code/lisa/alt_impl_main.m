%% setup; copied from state fdbk example
clear; close all; clc; 

% specify system matrices
sys    = LTISystem;
sys.Nx = 10;

alpha = 0.5; rho = 1; actDens = 0.8;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)]; % used in H2/HInf ctrl
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams       = SLSParams;
slsParams.tFIR_ = 20;
slsParams.obj_  = Objective.H2; % objective function

% simulation parameters
simParams           = SimParams;
simParams.tSim_     = 40;
simParams.w_        = zeros(sys.Nx, simParams.tSim_); % disturbance
simParams.w_(floor(sys.Nx/2), 1) = 10;
simParams.openLoop_ = false;
slsParams.mode_     = SLSMode.Basic;

slsOuts = state_fdbk_sls(sys, slsParams);
disp('Done')

%% find new implementations Rc, Mc; calculate stats
Tc          = 3;
slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc);
TMax        = max(Tc, slsParams.tFIR_);

% calculate stats
R  = slsOuts.R_; Rc = slsOuts_alt.R_;
M  = slsOuts.M_; Mc = slsOuts_alt.M_;

% make Rc, R equal length for easier comparisons
for t=slsParams.tFIR_+1:TMax % pad R, M with zeros if needed
    R{t} = zeros(sys.Nx, sys.Nx);
    M{t} = zeros(sys.Nu, sys.Nx);
end
for t=Tc+1:TMax % pad Rc, Mc with zeros if needed
    Rc{t} = zeros(sys.Nx, sys.Nx);
    Mc{t} = zeros(sys.Nu, sys.Nx);
end

Rdiff = 0; Mdiff = 0;
for t=1:TMax
    Rdiff = Rdiff + norm(full(R{t} - Rc{t}));
    Mdiff = Mdiff + norm(full(M{t} - Mc{t}));
end

Rdiff
Mdiff

slsParams_alt       = copy(slsParams);
slsParams_alt.tFIR_ = Tc;

[xOld, uOld] = simulate_system(sys, slsParams, slsOuts, simParams);
[xNew, uNew] = simulate_system(sys, slsParams_alt, slsOuts_alt, simParams);

xDiff = norm(xNew-xOld)
uDiff = norm(uNew-uOld)

plot_heat_map(xOld, sys.B2*uOld, 'Old');
plot_heat_map(xNew, sys.B2*uNew, 'New');
