classdef OptimizedCtrller < matlab.mixin.Copyable
    % Contains outputs (Rc, Mc) for 
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
      function obj = OptimizedCtrller()
        % initialize to zero instead of empty array
        obj.Rc_ = 0; obj.Mc_ = 0;
      end
    end
    
end