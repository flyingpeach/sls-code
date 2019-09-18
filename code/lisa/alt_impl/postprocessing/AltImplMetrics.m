classdef AltImplMetrics < matlab.mixin.Copyable
    % Contains metrics for alternative implementations Rc, Mc
    % Inherits handle class with deep copy functionality
    
    % R/M are original CL responses
    % Rc/Mc are new CL implementations
    % RTc/MTc are original CL responses clipped to be Tc long (if Tc<T)

    properties
      Tcs;              % list of Tcs (if sweeping over Tcs) or a single Tc
      sweepParams;      % if not sweeping over Tcs, the other things we sweep over
      sweepParamName;   % name of the sweep parameter (i.e. 'Tc')   
            
      tol; % tolerance used for zero and for rank calculations

      % properties of original R, M
      L1NormOrig;         % L1 norm of R/M
      RNonzero; MNonzero; % # nonzero entries (> zThresh) of R, M
      
      % properties of Rc, Mc
      L1Norms;                % L1 norms of Rc/Mc
      RcNonzeros; McNonzeros; % # nonzero entries (> zThresh) of Rc, Mc
      RDiffs;     MDiffs;     % sum over max(Tc, T) of ||Rc-R||, ||Mc-M||
      GcDiffs;    HcDiffs;    % Gc/Hc are new CL maps (using Rc, Mc) from w to x/u
                              % sum over max(Tc, T) of ||Gc-R||, ||Hc-M||
      
      % properties of RTc, MTc (identical to R, M when Tc >= T)
      L1NormsTc;                % L1 norms of RTc/MTc
      RTcNonzeros; MTcNonzeros; % # nonzero entries (> zThresh) of RTc, MTc
      RTcDiffs;    MTcDiffs;    % sum over max(Tc, T) of ||RTc-R||, ||MTc-M||
      GTcDiffs;    HTcDiffs;    % GTc/HTc are CL maps (using RTc, MTc) from w to x/u
                                % sum over max(Tc, T) of ||GTc-R||, ||HTc-M||
    end
    
    methods
      function obj = AltImplMetrics(tol, Tcs, sweepParamName, sweepParams) 
        % initializer
        if nargin == 2
            obj.sweepParamName = 'Tc';
        else
            obj.sweepParamName = sweepParamName;
            obj.sweepParams    = sweepParams;
        end
        
        obj.tol = tol;        
        obj.Tcs            = Tcs;
        
        if strcmp(obj.sweepParamName, 'Tc')
            numItems = length(Tcs);
        else
            numItems = length(sweepParams);
        end
        
        obj.L1NormOrig = 0;
        obj.RNonzero   = 0; obj.MNonzero = 0;
        
        obj.L1Norms    = zeros(numItems,1);
        obj.RcNonzeros = zeros(numItems,1); obj.McNonzeros = zeros(numItems,1);
        obj.RDiffs     = zeros(numItems,1); obj.MDiffs     = zeros(numItems,1);
        obj.GcDiffs    = zeros(numItems,1); obj.HcDiffs    = zeros(numItems,1);
        
        obj.L1NormsTc   = zeros(numItems,1);
        obj.RTcNonzeros = zeros(numItems,1); obj.MTcNonzeros = zeros(numItems,1);
        obj.RTcDiffs    = zeros(numItems,1); obj.MTcDiffs    = zeros(numItems,1);
        obj.GTcDiffs    = zeros(numItems,1); obj.HTcDiffs    = zeros(numItems,1); 
      end
    end
    
end