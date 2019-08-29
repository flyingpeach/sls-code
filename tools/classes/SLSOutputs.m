classdef SLSOutputs < matlab.mixin.Copyable
    % inherits handle class with deep copy functionality
    % contains outputs of SLS 
    % depending on the solver / mode called, not all outputs are set 
    
    properties 
      % based on (5.1),(5.2),(5.3)
      % reminder: the two disturbances are dx = B1*w, dy = D21*w
      
      R_; % dx to x transfer matrix (phi_xx)
      M_; % dx to u transfer matrix (phi_ux)

      % output feedback only
      N_; % dy to x transfer matrix (phi_xy)
      L_; % dy to u transfer matrix (phi_uy)
        
      clnorm_; % final (optimal value) of original objectiv
               % doesn't include regularization terms
      acts_;   % indices of actuators (u) kept after rfd
      
      % TODO: this enforces small gain on l1->l1; should generalize
      robustStab_; % inf norm of delta from (2.24), (4.22)
                   % <1 means we can guarantee stab
    end    
end