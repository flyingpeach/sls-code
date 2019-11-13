classdef SLSMode
    % Enumeration class containing options for sls closed-loop map 
    % design OR ctrller design (or joint design)

    enumeration     
      % Applicable to closed-loop map design OR ctrller design step      
      Basic % in ctrller design step, = minimize L1 norm of Rc/Mc
      Delayed
      Localized
      DAndL
      
      % Applicable only to ctrller design step
      % (Theoretically can use for closed-loop map but it's not currently
      %  implemented in that solver)
      EncourageDelay
      EncourageLocal
    end
end