function calc_lqr_costs(slsOutsCent, G, H)
% Outputs lqr cost of inf horizon, centralized, distributed, 
% and the reimplementation
% G, H are the CL maps (from w to x and w to u, respectively) 
% of the new implementation

%% Inf-horizon centralized
[P,L,K] = dare(full(sys.A), full(sys.B2), eye(sys.Nx), eye(sys.Nu));

infH2Cost = 0;
for i=1:sys.Nx % H2 cost is sum of all costs from init condn
    x0    = zeros(sys.Nx, 1);
    x0(i) = 1;
    infH2Cost = infH2Cost + x0'*P*x0;
end
infH2Cost

%% FIR centralized SLS
R = slsOuts.R_;
M = slsOuts.M_;
centH2Cost = 0;
for t=1:slsParams.tFIR_
    centH2Cost = centH2Cost + norm(full([R{t}; M{t}]), 'fro').^2;
end
centH2Cost

%% Distributed SLS
R = slsOuts.R_;
M = slsOuts.M_;
distH2Cost = 0;
for t=1:slsParams.tFIR_
    distH2Cost = distH2Cost + norm(full([R{t}; M{t}]), 'fro').^2;
end
distH2Cost

%% Reimplementation
tTotal = slsParams.tFIR_+50;
[G, H] = calc_cl_map(sys, slsParams, slsOuts, simParams, tTotal);
reImplH2Cost = 0;
for t=1:tTotal
    reImplH2Cost = reImplH2Cost + norm(full([G{t}; H{t}]), 'fro').^2;
end
reImplH2Cost

%% Output relative to infinite horizon
[sprintf('Centralized: %.3f', centH2Cost / infH2Cost);
 sprintf('Distributed: %.3f', distH2Cost / infH2Cost);
 sprintf('Reimplement: %.3f', reImplH2Cost / infH2Cost);]
