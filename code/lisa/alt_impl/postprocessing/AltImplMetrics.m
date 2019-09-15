classdef AltImplMetrics < matlab.mixin.Copyable
    % Contains metrics for alternative implementations Rc, Mc
    % Inherits handle class with deep copy functionality
    
    % R/M are original CL responses
    % Rc/Mc are new CL implementations
    % RTc/MTc are original CL responses clipped to be Tc long (if Tc<T)

    properties
      Tcs; % list of Tcs for which Rc/Mc was generated
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
      function obj = AltImplMetrics(Tcs, tol) % initializer
        obj.Tcs = Tcs;
        obj.tol = tol;
        numTcs  = length(Tcs);
        
        obj.L1NormOrig = 0;
        obj.RNonzero   = 0; obj.MNonzero = 0;
        
        obj.L1Norms    = zeros(numTcs,1);
        obj.RcNonzeros = zeros(numTcs,1); obj.McNonzeros = zeros(numTcs,1);
        obj.RDiffs     = zeros(numTcs,1); obj.MDiffs     = zeros(numTcs,1);
        obj.GcDiffs    = zeros(numTcs,1); obj.HcDiffs    = zeros(numTcs,1);
        
        obj.L1NormsTc   = zeros(numTcs,1);
        obj.RTcNonzeros = zeros(numTcs,1); obj.MTcNonzeros = zeros(numTcs,1);
        obj.RTcDiffs    = zeros(numTcs,1); obj.MTcDiffs    = zeros(numTcs,1);
        obj.GTcDiffs    = zeros(numTcs,1); obj.HTcDiffs    = zeros(numTcs,1); 
      end
    end
    
end