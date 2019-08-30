classdef SLSParams < matlab.mixin.Copyable
    % inherits handle class with deep copy functionality
    % contains parameters for SLS 
    % note that depending on the solver called, not all params may be used
    
    properties  
      mode_;     % an SLSMode (i.e. Basic, DLocalized)
      rfd_;      % whether we want rfd
      
      % used for all modes
      tFIR_;     % finite impulse response horizon
      
      % used for DLocalized and ApproxDLocalized modes only
      actDelay_; % actuation delay
      cSpeed_;   % communication speed
      d_;        % d-hop locality constraint
      
      robCoeff_; % regularization coeff for robust stability (used in approx-d sls only)
      rfdCoeff_; % regularization coeff for rfd (used in rfd only)
      
      obj_;      % an Objective (i.e. H2, HInf)
      xDes_;     % desired trajectory (only if obj is TrajTrack)
                 % do not directly set; use method setDesiredTraj
    end
    
    methods
      function obj = SLSParams()
        % initialize to zero instead of empty array
        obj.mode_     = 0; 
        obj.rfd_      = 0;
        
        obj.tFIR_     = 0;        
        obj.actDelay_ = 0;
        obj.cSpeed_   = 0;
        obj.d_        = 0;
        
        obj.robCoeff_ = 0; 
        obj.rfdCoeff_ = 0;
        
        obj.obj_      = 0;
        obj.xDes_     = 0;        
      end
        
      function setDesiredTraj(obj, xDes)
        % TODO: this is hacky
        % pad first column with zeros
        % also pad with zeros if xDes not specified for all time
        timediff  = obj.tFIR_ - size(xDes, 2);
        obj.xDes_ = [zeros(size(xDes, 1), 1), xDes, zeros(size(xDes, 1), timediff-1)];
      end
    end     
end