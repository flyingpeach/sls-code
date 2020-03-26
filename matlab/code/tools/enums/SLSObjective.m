classdef SLSObjective
    % Enumeration class containing options for objectives
    
    enumeration
      % Standard control objectives
      H2
      HInf
      L1
      
      % Regularizers
      Stability % for robust SLS
      RFD       % actuator regularization
      Locality  % penalize nonlocal patterns
      Delay     % penalize patterns that require fast communication
      
    end
end