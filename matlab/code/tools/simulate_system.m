function [x, u] = simulate_system(sys, simParams, Rc, Mc)
% Simulate system as per equation (2.8)
% Returns 
%    x, u      : state and actuation values
% Inputs
%    sys       : LTISystem containing system matrices
%    simParams : SimParams; parameters for the simulation
%    Rc, Mc    : controller implementation matrices
%                (can directly use R, M from SLSOutputs)

simParams.sanity_check();
        
x     = zeros(sys.Nx, simParams.tSim_); 
u     = zeros(sys.Nu, simParams.tSim_);
x_hat = zeros(sys.Nx, simParams.tSim_); 
w_hat = zeros(sys.Nx, simParams.tSim_);

T = length(Rc);

for t=1:1:simParams.tSim_-1
    if (simParams.openLoop_ ~= 1) % closed loop simulation
        for tau=1:1:min(t-1, T)
           u(:,t) = u(:,t) + Mc{tau}*w_hat(:,t-tau);
        end

        for tau=1:1:min(t-1, T-1)
           x_hat(:,t+1) = x_hat(:,t+1) + Rc{tau+1}*w_hat(:,t-tau);       
        end 
    end
    
    x(:,t+1) = sys.A*x(:,t) + sys.B1*simParams.w_(:,t)+ sys.B2*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end
