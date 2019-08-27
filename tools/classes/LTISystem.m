classdef LTISystem<handle
    % contains all matrices of an LTI system as per (3.1)
    properties
      % matrices
      A;  B1;  B2;  % state       : x(t+1)= A*x(t)  + B1*w(t)  + B2*u(t)
      C1; D11; D12; % reg output  : z_(t) = C1*x(t) + D11*w(t) + D12*u(t)
      C2; D21; D22; % measurement : y(t)  = C2*x(t) + D21*w(t) + D22*u(t)
      
      Nx; Nu; % number of states and number of actuators
    end
end