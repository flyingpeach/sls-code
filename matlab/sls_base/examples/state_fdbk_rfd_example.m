% note: run sequentially; some params from later sections depend on
% params specified in previous sections
clear; close all; clc;

% generate sys.A, sys.B2
Nx = 10; rho = 0.8; actDens = 1; randn('seed', 0);
sys = generate_rand_chain(Nx, rho, actDens);

sys.Nw  = sys.Nx; 
sys.Nz  = sys.Nu + sys.Nx;
sys.B1  = eye(sys.Nx); 
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D11 = sparse(sys.Nz, sys.Nw);
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];
sys.sanity_check();

% sls parameters (without rfd)
slsParams      = SLSParams();
slsParams.T_   = 15;
slsParams.add_objective(SLSObjective.H2, 1);
slsParams.add_constraint(SLSConstraint.ActDelay, 1);
slsParams.add_constraint(SLSConstraint.CommSpeed, 3);
slsParams.add_constraint(SLSConstraint.Locality, 4);

% sls parameters (with rfd)
slsParamsRFD = copy(slsParams);

% delayed and localized sls with rfd (sweep over different rfd coeffs)
powers    = -2:1:3;
numPowers = length(powers);

numActuators = zeros(numPowers, 1);
h2Objectives = zeros(numPowers, 1);

for i=1:numPowers
    slsParamsRFD.add_objective(SLSObjective.RFD, 10^powers(i));
    clMaps = state_fdbk_sls(sys, slsParamsRFD);
    
    % update system actuation
    sysAfterRFD = update_actuation(sys, clMaps);
    clMapsAfter = state_fdbk_sls(sysAfterRFD, slsParams);

    numActuators(i) = sysAfterRFD.Nu;
    objMtx          = get_objective_mtx(sysAfterRFD, clMapsAfter.R_, clMapsAfter.M_);
    h2Objectives(i) = get_H2_norm(objMtx);    
end

figure;
plot(numActuators, h2Objectives,'*-');
xlabel('# actuators'); ylabel('H2 objective val');
title('delayed + localized RFD tradeoff curve');
