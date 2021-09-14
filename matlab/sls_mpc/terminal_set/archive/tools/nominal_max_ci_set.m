function [Hx, hx, iters] = nominal_max_ci_set(sys, params)
% Gives terminal set (i.e. max. control invariant set)
% for some set defined by state and input constraints in params
% Output set is defined as x: Hx*x <= hx
% iters is the number of iterations required to compute the terminal set

maxIters = 50;
ZAB = get_sls_constraint(sys, 2);

% state-only constraints
Hx = [params.stateConsMtx_; -params.stateConsMtx_];
hx = [params.stateUB_; -params.stateLB_];

% input-only constraints
Hu = [params.inputConsMtx_; -params.inputConsMtx_];
hu = [params.inputUB_; -params.inputLB_];

% subsequent steps are Pre(Omega_k) \intersect (Omega_k)
for iters=1:maxIters
    iters
     last_Hx = Hx;
     last_hx = hx;
               
     % Hx, hx are for state only, but to generate the next set we need 
     % to include input constraints as well since we will be doing
     % H*Phi*x which includes Phix(1), Phix(2), and Phiu(1)     
     Hxu = blkdiag(Hx, Hx, Hu);
     hxu = [hx; hx; hu];
     
     if isempty(Hu) % Fix dimensions
        Hxu = [Hxu zeros(size(Hxu, 1), sys.Nu)];
     end
     
     [Hx, hx] = nominal_set(Hxu, ZAB, hxu, sys.Nx);
  
     % Check for convergence
     lastP = Polyhedron(last_Hx, last_hx); 
     thisP = Polyhedron(Hx, hx);
     if lastP == thisP
        break
     end
end

thisP.minHRep();
Hx = thisP.A; hx = thisP.b;

if iters == maxIters    
    fprintf('Terminal set computation reached %d iters and did not converge\n', maxIters);
end

end