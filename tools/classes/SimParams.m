classdef SimParams < matlab.mixin.Copyable
    % inherits handle class with deep copy functionality
    % contains parameters for simulating the system

    properties  
      tSim_;     % amount of time to simulate
      w_;        % disturbance w(t) (as per (3.1)
      openLoop_; % whether to simulate the system open loop only
    end
    
    methods
      function obj = SimParams()
        % initialize to zero instead of empty array
        obj.tSim_     = 0; 
        obj.w_        = 0; 
        obj.openLoop_ = 0;
      end
        
      function txt = print_and_check(obj)
        modestr = 'closed-loop';
        if obj.openLoop_ == 1
            modestr = 'open-loop';
        end
        txt = ['tSim=', num2str(obj.tSim_), ', ', modestr];
        
        % sanity check
        if size(obj.w_, 2) < obj.tSim_
            error('[SLS ERROR] The specified length of the disturbance (w) is less than tSim!')
        end
      end
    end
end