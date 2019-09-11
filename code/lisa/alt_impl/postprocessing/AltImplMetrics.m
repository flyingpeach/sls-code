classdef AltImplMetrics < matlab.mixin.Copyable
    % Contains metrics for alternative implementations Rc, Mc
    % Inherits handle class with deep copy functionality
    
    % R/M are original CL responses
    % Rc/Mc are new CL responses
    % RTc/MTc are original CL responses clipped to be Tc long (if Tc<T)

    properties
      Tcs;     % list of Tcs for which Rc/Mc was generated
      zThresh; % anything lower than this is counted as 0 for RNonzero, etc.

      % properties of original R, M
      L1NormOrig;         % L1 norm of R/M
      RNonzero; MNonzero; % # nonzero entries (> zThresh) of R, M
      
      % properties of Rc, Mc
      L1Norms;                % L1 norms of Rc/Mc
      RcNonzeros; McNonzeros; % # nonzero entries (> zThresh) of Rc, Mc
      RDiffs;     MDiffs;     % sum over max(Tc, T) of ||Rc-R||, ||Mc-M||
      xDiffs;     uDiffs;     % differences btw CL x, u of Rc/Mc and R/M
      
      % properties of RTc, MTc (identical to R, M when Tc >= T)
      L1NormsTc;                % L1 norms of RTc/MTc
      RTcNonzeros; MTcNonzeros; % # nonzero entries (> zThresh) of RTc, MTc
      RTcDiffs; MTcDiffs;       % sum over max(Tc, T) of ||RTc-R||, ||MTc-M||
      xTcDiffs; uTcDiffs;       % differences btw CL x, u of RTc/MTc and R/M  
    end
    
    methods
      function obj = AltImplMetrics(Tcs, zThresh) % initializer
        obj.Tcs     = Tcs;
        obj.zThresh = zThresh;
        numTcs      = length(Tcs);
        
        obj.L1NormOrig = 0;
        obj.RNonzero   = 0; obj.MNonzero = 0;
        
        obj.L1Norms    = zeros(numTcs,1);
        obj.RcNonzeros = zeros(numTcs,1); obj.McNonzeros = zeros(numTcs,1);
        obj.RDiffs     = zeros(numTcs,1); obj.MDiffs     = zeros(numTcs,1);
        obj.xDiffs     = zeros(numTcs,1); obj.uDiffs     = zeros(numTcs,1);
        
        obj.L1NormsTc   = zeros(numTcs,1);
        obj.RTcNonzeros = zeros(numTcs,1); obj.MTcNonzeros = zeros(numTcs,1);
        obj.RTcDiffs    = zeros(numTcs,1); obj.MTcDiffs    = zeros(numTcs,1);
        obj.xTcDiffs    = zeros(numTcs,1); obj.uTcDiffs    = zeros(numTcs,1);              
      end
    end
    
end