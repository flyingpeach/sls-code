function [x, u] = simulate_system(sys, slsParams, slsOuts, TMax, w, openLoop)
% Simulate system as per equation (2.8)
% Returns 
%    x, u     : state and actuation values
% Inputs
%    sys       : LTISystem containing system matrices
%    slsParams : SLSParams containing parameters
%    slsOuts   : SLSOutputs containing system responses and other info
%    TMax      : total amount of time to simulate for
%    w         : disturbance values (Nx by TMax)
%    openLoop  : optional, set to 1 if you want open loop simulation

if nargin == 5
    openLoop = 0;
end
   
x     = zeros(sys.Nx, TMax); 
u     = zeros(sys.Nu, TMax);
x_hat = zeros(sys.Nx, TMax); 
w_hat = zeros(sys.Nx, TMax);

for t=1:1:TMax-1
    if (openLoop ~= 1) % closed loop simulation
        for tau=1:1:min(t-1, slsParams.tFIR_)
           u(:,t) = u(:,t) + slsOuts.M_{tau}*w_hat(:,t-tau);
        end

        for tau=1:1:min(t-1, slsParams.tFIR_-1)
           x_hat(:,t+1) = x_hat(:,t+1) + slsOuts.R_{tau+1}*w_hat(:,t-tau);       
        end 
    end
    
    x(:,t+1) = sys.A*x(:,t) + sys.B1*w(:,t)+ sys.B2*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end
