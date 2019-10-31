classdef SLSMode
    % Enumeration class containing modes of sls
    enumeration
      Basic, 
      Localized,
      Delayed,
      DAndL, % delayed and localized
      ApproxDAndL % approximately delayed and localized
    end
end