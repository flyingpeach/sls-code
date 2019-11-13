classdef Ctrller < matlab.mixin.Copyable
    % Contains controller implementation matrices (Rc, Mc)
    % Inherits handle class with deep copy functionality
    
    properties
      Rc_;
      Mc_;
    end
    
    methods
      function obj = Ctrller()
        % initialize to zero instead of empty array
        obj.Rc_ = 0; obj.Mc_ = 0;
      end
    end
    
end