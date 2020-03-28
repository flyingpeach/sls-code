function [G, H] = get_cl_map(sys, ctrller, tTotal)
simParams           = SimParams();
simParams.tSim_     = tTotal + 1; % reactions to disturbance will start @ t=2
simParams.openLoop_ = false;

for t=1:tTotal
    G{t} = zeros(sys.Nx, sys.Nx); % CL map from w to x
    H{t} = zeros(sys.Nu, sys.Nx); % CL map from w to u
end

for i=1:sys.Nx
    simParams.w_ = zeros(sys.Nx, simParams.tSim_);
    simParams.w_(i, 1) = 1;

    [x, u]  = simulate_system(sys, simParams, ctrller.Rc_, ctrller.Mc_);

    for t=1:tTotal
        G{t}(:,i) = x(:,t+1);
        H{t}(:,i) = u(:,t+1);
    end
end
