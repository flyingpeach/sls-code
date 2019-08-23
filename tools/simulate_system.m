function [x, u] = simulate_system(sys, R, M, w, TFIR, TMax)
% Simulate system as per equation (2.8)
% Returns 
%    x, u : state and actuation values
% Inputs
%    sys  : an LTISystem
%    R, M : system response designed by SLS
%    w    : disturbance values (Nx by TMax)  
%    TFIR : finite impulse response time
%    TMax : total amount of time to simulate for

x     = zeros(sys.Nx, TMax); 
u     = zeros(sys.Nu, TMax);
x_hat = zeros(sys.Nx, TMax); 
w_hat = zeros(sys.Nx, TMax);

for t=1:1:TMax-1
    for tau=1:1:min(t-1, TFIR)
       u(:,t) = u(:,t) + M{tau}*w_hat(:,t-tau);
    end

    for tau=1:1:min(t-1, TFIR-1)
       x_hat(:,t+1) = x_hat(:,t+1) + R{tau+1}*w_hat(:,t-tau);       
    end    

    x(:,t+1) = sys.A*x(:,t) + sys.B1*w(:,t)+ sys.B2*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end
