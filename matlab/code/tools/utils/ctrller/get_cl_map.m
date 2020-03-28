function [R_, M_] = get_cl_map(sys, ctrller, tTotal)
% Obtains implemented closed-loop maps (R_, M_) via simulation
% Results will be inaccurate if ctrller unstable or tTotal too small

simParams       = SimParams();
simParams.tSim_ = tTotal + 1; % reactions to disturbance will start @ t=2

R_ = cell(tTotal, 1); 
M_ = cell(tTotal, 1);

for t=1:tTotal
    R_{t} = zeros(sys.Nx, sys.Nx); % CL map from w to x
    M_{t} = zeros(sys.Nu, sys.Nx); % CL map from w to u
end

for i=1:sys.Nx
    simParams.w_ = zeros(sys.Nx, simParams.tSim_);
    simParams.w_(i, 1) = 1;

    [x, u]  = simulate_system(sys, ctrller, simParams);

    for t=1:tTotal
        R_{t}(:,i) = x(:,t+1);
        M_{t}(:,i) = u(:,t+1);
    end
end
