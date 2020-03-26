function [x, u] = simulate_state_fdbk(sys, ctrller, params)
% Simulate state feedback system as per equation (2.8)
% Returns 
%    x, u    : state and actuation values
% Inputs
%    sys     : LTISystem containing system matrices
%    ctrller : Ctrller containing implementation matrices
%    params  : SimParams containing parameters for simulation
params.sanity_check();

tSim = params.tSim_;
Rc   = ctrller.Rc_;
Mc   = ctrller.Mc_;

x     = zeros(sys.Nx, tSim); 
u     = zeros(sys.Nu, tSim);
x_hat = zeros(sys.Nx, tSim);
w_hat = zeros(sys.Nx, tSim);

T = length(Rc);

for t=1:1:tSim-1
    for tau=1:1:min(t-1, T)
       u(:,t) = u(:,t) + Mc{tau}*w_hat(:,t-tau);
    end

    for tau=1:1:min(t-1, T-1)
       x_hat(:,t+1) = x_hat(:,t+1) + Rc{tau+1}*w_hat(:,t-tau);       
    end 
    
    x(:,t+1) = sys.A*x(:,t) + sys.B1*params.w_(:,t)+ sys.B2*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end
