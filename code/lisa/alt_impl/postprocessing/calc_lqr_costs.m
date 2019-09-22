function calc_lqr_costs(slsOutsCent, G, H)
% Outputs lqr cost of inf horizon, centralized, distributed, 
% and the reimplementation
% G, H are the CL maps (from w to x and w to u, respectively) 
% of the new implementation


% Inf-horizon centralized
[P,L,K] = dare(full(sys.A), full(sys.B2), eye(sys.Nx), eye(sys.Nu));

infH2Cost = 0;
for i=1:sys.Nx % H2 cost is sum of all costs from init condn
    x0    = zeros(sys.Nx, 1);
    x0(i) = 1;
    infH2Cost = infH2Cost + x0'*P*x0;
end

% FIR centralized SLS
centH2Cost = 0;
RCent = slsOutsCent.R_;
MCent = slsOutsCent.M_;
for t=1:slsParams.tFIR_
    centH2Cost = centH2Cost + norm(full([RCent{t}; MCent{t}]), 'fro').^2;
end

% Distributed SLS (for delay=1, setup1) 
distH2Cost = 16.7383;

% Reimplementation
for t=1:slsParams.tFIR_
    reImplH2Cost = reImplH2Cost + norm(full([G{t}; H{t}]), 'fro').^2;
end

% Output relative to infinite horizon
[sprintf('Centralized: %.4f', centH2Cost / infH2Cost);
 sprintf('Distributed: %.4f', distH2Cost / infH2Cost);
 sprintf('Reimplement: %.4f', reImplH2Cost / infH2Cost);]
