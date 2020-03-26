classdef LTISystem < matlab.mixin.Copyable
    % Contains all matrices of an LTI system as per (3.1)
    % Inherits handle class with deep copy functionality

    properties
      % matrices
      A;  B1;  B2;  % state       : x(t+1)= A*x(t)  + B1*w(t)  + B2*u(t)
      C1; D11; D12; % reg output  : z_(t) = C1*x(t) + D11*w(t) + D12*u(t)
      C2; D21; D22; % measurement : y(t)  = C2*x(t) + D21*w(t) + D22*u(t)
      
      Nx; % # states
      Nw; % # disturbances
      Nu; % # actuators
      Nz; % # regulated outputs
      Ny; % # outputs      
    end
    
    methods
      function obj = LTISystem()
        obj.A  = []; obj.B1  = []; obj.B2  = []; 
        obj.C1 = []; obj.D11 = []; obj.D12 = [];
        obj.C2 = []; obj.D21 = []; obj.D22 = [];
        
        obj.Nx = 0; 
        obj.Nw = 0;
        obj.Nu = 0;
        obj.Nz = 0;
        obj.Ny = 0;
      end

      function obj = updateActuation(oldObj, slsOuts)
        % make a new system with the dynamics of the old system and updated
        % actuation (based on rfd output)
        obj     = copy(oldObj);
        obj.B2  = oldObj.B2(:, slsOuts.acts_);
        obj.D12 = oldObj.D12(:, slsOuts.acts_);
        obj.Nu  = size(slsOuts.acts_, 1);      
      end
      
      function sanityCheck(obj)
        % sanity check that dimensions are consistent
        if size(obj.A,1) ~= obj.Nx || size(obj.A,2) ~= obj.Nx
            sls_error('A matrix has incorrect dimensions')
        end
        if size(obj.B1,1) ~= obj.Nx || size(obj.B1,2) ~= obj.Nw
            sls_error('B1 matrix has incorrect dimensions')
        end
        if size(obj.B2,1) ~= obj.Nx || size(obj.B2,2) ~= obj.Nu
            sls_error('B2 matrix has incorrect dimensions')
        end

        if size(obj.C1,1) ~= obj.Nz || size(obj.C1,2) ~= obj.Nx
            sls_error('C1 matrix has incorrect dimensions')
        end
        if size(obj.D11,1) ~= obj.Nz || size(obj.D11,2) ~= obj.Nw
            sls_error('D11 matrix has incorrect dimensions')
        end
        if size(obj.D12,1) ~= obj.Nz || size(obj.D12,2) ~= obj.Nu
            sls_error('D12 matrix has incorrect dimensions')
        end

        if size(obj.C2,1) ~= obj.Ny || size(obj.C2,2) ~= obj.Nx
            sls_error('C2 matrix has incorrect dimensions')
        end
        if size(obj.D21,1) ~= obj.Ny || size(obj.D21,2) ~= obj.Nw
            sls_error('D21 matrix has incorrect dimensions')
        end
        if size(obj.D22,1) ~= obj.Ny || size(obj.D22,2) ~= obj.Nu
            sls_error('D22 matrix has incorrect dimensions')
        end
    end

end