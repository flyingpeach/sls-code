classdef SLSOutputs < matlab.mixin.Copyable
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
        
      clnorm_; % final (optimal value) of original objective
               % doesn't include regularization terms
      acts_;   % indices of actuators (u) kept after rfd
      
      solveStatus_; % cvx_status
      
      % TODO: this enforces small gain on l1->l1; should generalize
      robustStab_; % inf norm of delta from (2.24), (4.22)
                   % <1 means we can guarantee stab
    end
    
    methods
      function obj = SLSOutputs()
        % initialize to zero instead of empty array
        obj.R_ = 0; obj.M_ = 0; obj.N_ = 0; obj.L_ = 0;
        obj.clnorm_     = 0;
        obj.acts_       = 0;
        obj.robustStab_ = 0;
      end
    end
    
end