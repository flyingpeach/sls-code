classdef CtrllerStats < matlab.mixin.Copyable
    % Contains stats for controller matrices Rc, Mc
    % Inherits handle class with deep copy functionality
    
    % R/M are original CL responses
    % Rc/Mc are new CL implementations
    % RTc/MTc are original CL responses clipped to be Tc long (if Tc<T)

    properties
      sweepParams;      % list of parameter values (optional)
      sweepParamName;   % name of the sweep parameter (i.e. 'delay')   
            
      % properties of original R, M
      L1NormOrig;         % L1 norm of R/M
      IntSpecRadiusOrig;  % original spectral radius of internal dynamics
      LQRCostOrig;        % original LQR cost
      
      % properties of Rc, Mc
      L1Norms;                % L1 norms of Rc/Mc
      GcDiffs;    HcDiffs;    % Gc/Hc are new CL maps (using Rc, Mc) from w to x/u
                              % sum over max(Tc, T) of ||Gc-R||, ||Hc-M||
      IntSpecRadii_c;         % spectral radii of internal dynamics (< 1 means stable)
      LQRCosts;               % LQR costs of new implementations
    end
    
    methods
      function obj = CtrllerStats(sweepParamName, sweepParams) 
        % initializer
        if nargin == 3
            obj.sweepParamName = sweepParamName;
            obj.sweepParams    = sweepParams;
            numItems           = length(sweepParams);
        else
            obj.sweepParamName = 0;
            obj.sweepParams    = 0;
            numItems           = 1;
        end

        obj.L1NormOrig = 0;
        obj.IntSpecRadiusOrig = 0;
        obj.LQRCostOrig = 0;
        
        obj.L1Norms    = zeros(numItems,1);
        obj.GcDiffs    = zeros(numItems,1); obj.HcDiffs    = zeros(numItems,1);
        obj.IntSpecRadii_c = zeros(numItems,1);
        obj.LQRCosts = zeros(numItems,1);
      end
    end
    
end