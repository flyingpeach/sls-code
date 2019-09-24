classdef AltImplMode
    % Enumeration class containing alternative implementation modes
    enumeration
      ExactOpt % Minimizes L1 norm of [Rc; Mc]
               % by searching over Rc_specific + Rc_homogenous

               % note that this is implemented in an ApproxLS way to bypass
               % numerical issues

      Analytic % [Rc; Mc] = -F2 \ F1

      ApproxLS % ExplicitOpt with Rc_specific = LS solution and 
               % Rc_homogenous including right singular vectors 
               % corresponding to small nonzero singular values

      ApproxLeaky % ExplicitOpt with penalized leaky term
      
      StrictDelay % ExplicitOpt with strict comm delay constraints and M1 relaxation
      
      StrictLocal % ExplicitOpt with strict locality constraint and M1 relaxation
      
      EncourageDelay % ExplicitOpt with optimization encouraging tolerance
                     % for communication delays
    end
end