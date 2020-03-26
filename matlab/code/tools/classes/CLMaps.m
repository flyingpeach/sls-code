classdef CLMaps < matlab.mixin.Copyable
    % Contains outputs of SLS
    % Depending on the solver / mode called, not all outputs are set 
    % Inherits handle class with deep copy functionality
    
    properties 
      % based on (5.1),(5.2),(5.3)
      % reminder: the two disturbances are dx = B1*w, dy = D21*w
      
      R_; % dx to x transfer matrix (phi_xx)
      M_; % dx to u transfer matrix (phi_ux)

      % output feedback only
      N_; % dy to x transfer matrix (phi_xy)
      L_; % dy to u transfer matrix (phi_uy)

      acts_;        % indices of actuators (u) kept after rfd
      solveStatus_; % cvx_status
      
    end    
end