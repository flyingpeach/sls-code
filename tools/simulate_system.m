function [x, u] = simulate_system(A, B, B1, R, M, w, TFIR, TMax)
% Simulate system as per equation (2.8)
% Returns 
%    x, u     : state and actuation values
% Inputs
%    A, B, B1 : system matrices as per equation (3.1)
%    R, M     : system response designed by SLS
%    w        : disturbance values (Nx by TMax)  
%    TFIR     : finite impulse response time
%    TMax     : total amount of time to simulate for

Nx = size(A, 1); Nu = size(B, 2);

x     = zeros(Nx,TMax); u     = zeros(Nu,TMax);
x_hat = zeros(Nx,TMax); w_hat = zeros(Nx,TMax);

for t=1:1:TMax-1
    for tau=1:1:min(t-1, TFIR)
       u(:,t) = u(:,t) + M{tau}*w_hat(:,t-tau);
    end

    for tau=1:1:min(t-1, TFIR-1)
       x_hat(:,t+1) = x_hat(:,t+1) + R{tau+1}*w_hat(:,t-tau);       
    end    

    x(:,t+1) = A*x(:,t) + B1*w(:,t)+ B*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end
