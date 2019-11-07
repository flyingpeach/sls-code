classdef CtrllerMode
    % Enumeration class containing alternative implementation modes
    enumeration
      OptL1 % Minimizes L1 norm of [Rc; Mc] 
            % by searching over Rc_specific + Rc_homogenous
            % Rc_homogenous including right singular vectors 
            % corresponding to small nonzero singular values
     
      OptL1Delayed   % OptL1 with delay constraints only
      OptL1Localized % OptL1 with locality constraints only
      OptL1DAndL     % OptL1 with delay and locality constraints
     
      EncourageDelay % OptL1 that encourages tolerance for communication delay
      EncourageLocal % OptL1 that encourages tolerance for locality
    end
end