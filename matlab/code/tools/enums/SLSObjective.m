classdef SLSObjective
    % Enumeration class containing options for objectives
    
    enumeration
      % Standard control objectives
      H2
      HInf
      L1
      
      % Regularizers
      RFD       % actuator regularization
      
    end
end