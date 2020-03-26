classdef SimParams < matlab.mixin.Copyable
    % Contains parameters for simulating the system
    % Inherits handle class with deep copy functionality

    properties  
      tSim_;     % amount of time to simulate
      w_;        % disturbance w(t) (as per (3.1))
    end
    
    methods
      function obj = SimParams()
        % initialize to zero instead of empty array
        obj.tSim_     = 0; 
        obj.w_        = 0; 
      end
        
      function sanity_check(obj)
        if size(obj.w_, 2) < obj.tSim_
            sls_error('w_ must be at least as long as tSim')
        end
      end  
    end
    
end