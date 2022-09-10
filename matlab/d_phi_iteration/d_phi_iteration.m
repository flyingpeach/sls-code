function [Phi_xs, Phi_us, betas, nomCosts] = d_phi_iteration(sys, params)
% Implements D-Phi iteration for state feedback as presented in https://arxiv.org/abs/2204.02493
% Inputs
%   sys    : LTISystem containing relevant system matrices
%   params : DPhiParams class containing relevant parameters
% Outputs
%   Phi_xs, Phi_us : list of Phi_x and Phi_u (one for each iteration)
%   betas          : robust stability margin, one per iteration
%   nomCosts       : nominal performance cost, one per iteration

params.sanity_check()
Nx = size(sys.A, 1); MAX_ITERS = 100;

if ~params.randomizeD_
    if ismember(params.stabType_, [StabType.L1, StabType.LInf])
        fprintf('Note: this procedure is not separable in the D step\n');
    end
end

% Initialize
D      = eye(Nx);
Phi_xs = {0}; Phi_us = {0}; betas = [Inf]; nomCosts = [Inf];

fprintf('%s robust control\n', params.stabType_);
fprintf('\tStep size for robust bound: %.2e\n', params.betaStep_);
fprintf('\tStop when robust bound reaches: %.2e\n', params.betaStop_);
fprintf('Starting D-Phi iterations');
if params.randomizeD_
    fprintf(' with randomizing D step...\n');
else
    fprintf(' with minimizing D step...\n');
end

for k=2:MAX_ITERS
    betas(k) = betas(k-1) - params.betaStep_;
    [Phi_xs{k}, Phi_us{k}, nomCosts(k)] = phi_step(sys, params, D, betas(k));
    
    if isnan(nomCosts(k)) % phi step was infeasible
        exitFlag = true;
        if params.randomizeD_ % try to use D to advance beta instead
            D = d_step_randomize(M, betas(k), params.stabType_);
            if ~isnan(D(1,1))
                [Phi_xs{k}, Phi_us{k}, nomCosts(k)] = phi_step(sys, params, D, betas(k));
                M = make_m(params, Phi_xs{k}, Phi_us{k});
                betas(k) = get_bound(D, M, params.stabType_);
                exitFlag = false;
            end
        end
        if exitFlag
            % Last iteration was not feasible; don't return NaN values
            Phi_xs(end) = []; Phi_us(end) = []; betas(end) = []; nomCosts(end) = [];
            fprintf('\nPhi step infeasible. \nCannot improve robust bound anymore. Exiting.\n');
            break;       
        end
    end
    
    M = make_m(params, Phi_xs{k}, Phi_us{k});
    
    if params.randomizeD_ % Algorithm 2, randomizing D step
        betas(k) = get_bound(D, M, params.stabType_);
        D        = d_step_randomize(M, betas(k), params.stabType_);
    else % Algorithm 1, minimizing D step   
       [D, betas(k)] = d_step_minimize(M, params.stabType_);
    end
    
   fprintf('NomPerf: %+.2e | Bound: %+.2e\n', nomCosts(k), betas(k));
   if betas(k) <= params.betaStop_
       break; % reached desired robust bound, stop
   end
end

if k == MAX_ITERS
    error('\nD-Phi iteration limit of %d steps was reached. \nRobust bound may still be improvable. Exiting.\n', MAX_ITERS);
end

% First step is a placeholder; remove it
Phi_xs(1) = []; Phi_us(1) = []; betas(1) = []; nomCosts(1) = [];

end