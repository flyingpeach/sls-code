function [G, H] = calc_cl_map(sys, slsParams_alt, slsOuts_alt, simParams, tTotal)
simParams.tSim_ = tTotal + 1; % reactions to disturbance will start @ t=2

for t=1:tTotal
    G{t} = zeros(sys.Nx, sys.Nx); % CL map from w to x
    H{t} = zeros(sys.Nu, sys.Nx); % CL map from w to u
end

for i=1:sys.Nx
    simParams.w_ = zeros(sys.Nx, simParams.tSim_);
    simParams.w_(i, 1) = 1;
    [x, u]  = simulate_system(sys, slsParams_alt, slsOuts_alt, simParams);

    for t=1:tTotal
        G{t}(:,i) = x(:,t+1);
        H{t}(:,i) = u(:,t+1);
    end
end

