classdef SLSParams < matlab.mixin.Copyable
    % inherits handle class with deep copy functionality
    % contains parameters for SLS 
    % note that depending on the solver called, not all params may be used
    
    properties  
      actDelay_; % actuation delay
      cSpeed_;   % communication speed
      
      d_;        % d-hop locality constraint
      tFIR_;     % finite impulse response horizon
      
      robCoeff_; % regularization coeff for robust stability (used in approx-d sls only)
      rfdCoeff_; % regularization coeff for rfd (used in rfd only)
      
      mode_;     % an SLSMode (i.e. Basic, DLocalized)
      rfd_;      % whether we want rfd
      
      obj_;      % an Objective (i.e. H2, HInf)
      xDes_;     % desired trajectory (only if obj is TrajTrack)
                 % do not directly set; use method setDesiredTraj
    end
    
    methods
      function setDesiredTraj(obj, xDes)
        % TODO: pad first column with zeros (slightly hacky)
        % also pad with zeros if xDes not specified for all time
        timediff  = obj.tFIR_ - size(xDes, 2);
        obj.xDes_ = [zeros(size(xDes, 1), 1), xDes, zeros(size(xDes, 1), timediff-1)];
      end
    end     
end