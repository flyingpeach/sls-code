function [x, u] = simulate_system(sys, params, TMax, R, M, w, openLoop)
% Simulate system as per equation (2.8)
% Returns 
%    x, u     : state and actuation values
% Inputs
%    sys      : an LTISystem
%    params   : SLSParams containing parameters
%    TMax     : total amount of time to simulate for
%    R, M     : system response designed by SLS
%    w        : disturbance values (Nx by TMax)
%    openLoop : optional, set to 1 if you want open loop simulation

if nargin == 6
    openLoop = 0;
end
   
x     = zeros(sys.Nx, TMax); 
u     = zeros(sys.Nu, TMax);
x_hat = zeros(sys.Nx, TMax); 
w_hat = zeros(sys.Nx, TMax);

for t=1:1:TMax-1
    if (openLoop == 1) % open loop simulation
        u(:,t) = zeros(sys.Nu, 1);
    else % closed loop simulation
        for tau=1:1:min(t-1, params.tFIR_)
           u(:,t) = u(:,t) + M{tau}*w_hat(:,t-tau);
        end

        for tau=1:1:min(t-1, params.tFIR_-1)
           x_hat(:,t+1) = x_hat(:,t+1) + R{tau+1}*w_hat(:,t-tau);       
        end 
    end
    
    x(:,t+1) = sys.A*x(:,t) + sys.B1*w(:,t)+ sys.B2*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end
