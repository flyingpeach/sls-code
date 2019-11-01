classdef LTISystem < matlab.mixin.Copyable
    % Contains all matrices of an LTI system as per (3.1)
    % Inherits handle class with deep copy functionality

    properties
      % matrices
      A;  B1;  B2;  % state       : x(t+1)= A*x(t)  + B1*w(t)  + B2*u(t)
      C1; D11; D12; % reg output  : z_(t) = C1*x(t) + D11*w(t) + D12*u(t)
      C2; D21; D22; % measurement : y(t)  = C2*x(t) + D21*w(t) + D22*u(t)
      
      Nx; Nu; % number of states and number of actuators
    end
    
    methods
      function obj = LTISystem()
        % initialize to zero instead of empty array
        obj.A  = 0; obj.B1  = 0; obj.B2  = 0; 
        obj.C1 = 0; obj.D11 = 0; obj.D12 = 0;
        obj.C2 = 0; obj.D21 = 0; obj.D22 = 0;
        obj.Nx = 0; obj.Nu = 0;
      end

      function obj = updateActuation(oldObj, slsOuts)
        % make a new system with the dynamics of the old system and updated
        % actuation (based on rfd output)
        obj     = copy(oldObj);
        obj.B2  = oldObj.B2(:, slsOuts.acts_);
        obj.D12 = oldObj.D12(:, slsOuts.acts_);
        obj.Nu  = size(slsOuts.acts_, 1);      
      end
    end

end