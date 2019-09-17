classdef AltImplMode
    % Enumeration class containing alternative implementation modes
    enumeration
      ImplicitOpt % Minimizes L1 norm of [Rc; Mc] 
                  % with [Rc; Mc] == [R; M] * [zI-A -B] * [Rc; Mc] 

      ExplicitOpt % Minimizes L1 norm of [Rc; Mc]
                  % with F2 * [Rc; Mc] == -F1 
                  % where F1, F2 do not depend on [Rc, Mc]
                   
      Analytic    % [Rc; Mc] = -F2 \ F1

      ApproxDrop  % ExplicitOpt but with some dropped constraints

      ApproxLeaky % Minimizes L1 norm of [Rc; Mc] + weight * norm(Delta)
                  % where Delta = F2 * [Rc; Mc] + F1 (should be zero)
    end
end