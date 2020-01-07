%% solve original SLS problem
clear; close all; clc; 

% set up system / parameters
sys    = LTISystem();
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

slsParams       = SLSParams();
slsParams.T_    = 20;
slsParams.obj_  = Objective.H2;

% centralized SLS
slsParams.mode_ = SLSMode.Basic;
slsOuts         = state_fdbk_sls(sys, slsParams);

% infinite horizon
[P,L,K] = dare(full(sys.A), full(sys.B2), eye(sys.Nx), eye(sys.Nu));

infH2Cost = 0;
for i=1:sys.Nx % H2 cost is sum of all costs from init condition
    x0    = zeros(sys.Nx, 1);
    x0(i) = 1;
    infH2Cost = infH2Cost + x0'*P*x0;
end

%% sample space of implementation matrices
eps_base   = 2.22e-16;
eps_nullsp = eps_base.^(3/8);

cParams             = CtrllerParams();
cParams.mode_       = slsParams.mode_;
cParams.eps_nullsp_ = eps_nullsp;
cParams.T_          = slsParams.T_;

F  = get_ctrller_constraint(sys, slsOuts, cParams.T_);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

RMc_p         = F2 \ (-F1);
nullsp        = get_nullsp(F2, cParams.eps_nullsp_);
solnSpaceSize = size(nullsp, 2);

% draw from unif(-unifBounds, unifBounds)
unifBounds   = 1e3;  
numPts       = 1000;
l1norms      = zeros(numPts, 1);
specRads     = zeros(numPts, 1);
sanityChecks = zeros(numPts, 1); % check that it's still in the space

for i=1:numPts
    nullCoeff = (rand(solnSpaceSize, sys.Nx) - 0.5) * 2 * unifBounds;
    RMc = RMc_p + nullsp*nullCoeff;
    
    % expect a very small norm
    sanityChecks(i) = norm(F2 * RMc + F1);
    
    [Rc, Mc] = block_to_cell(RMc, cParams, sys);

    for t=1:cParams.T_ % calculate L1 norm
        l1norms(i) = l1norms(i) + norm([Rc{t}; Mc{t}], 1);
    end
        
    % calculate spec radius
    specRads(i) = check_int_stability(sys, Rc, Mc);
end

[Rc_p, Mc_p] = block_to_cell(RMc_p, cParams, sys);

l1normP = 0;
for t=1:cParams.T_
    l1normP = l1normP + norm([Rc_p{t}; Mc_p{t}], 1);
end

l1normP
specRadP = check_int_stability(sys, Rc_p, Mc_p)
norm(F2 * RMc_p + F1) % for comparison with sanity check value

figure(1)
subplot(2,1,1)
histogram(l1norms)
title('L1 norms')

subplot(2,1,2)
histogram(specRads)
title('Spec radii')

figure(2)
histogram(sanityChecks)
title('Sanity check (expect small)')

%% helper functions (copied from find_ctrller, local functions)

function range = get_range(idx, size)
% helper function to convert cells of block matrices into giant matrix
% copied from get_F (inline)
range = size*(idx-1)+1:size*(idx-1)+size;
end 

function [Rc, Mc] = block_to_cell(RMc, cParams, sys)
% Rc, Mc are cell structure (i.e. Rc{t})
% RMc is one stacked matrix (Rc first, then Mc)
Rcs = RMc(1:sys.Nx * (cParams.T_-1), :);
Mcs = RMc(sys.Nx * (cParams.T_-1) + 1:end, :);

for t = 1:cParams.T_
    tx    = get_range(t-1, sys.Nx);
    tu    = get_range(t, sys.Nu);
    Mc{t} = Mcs(tu,:);

    if t==1
        Rc{t} = eye(sys.Nx);
    else        
        Rc{t} = Rcs(tx,:);
    end
end
end