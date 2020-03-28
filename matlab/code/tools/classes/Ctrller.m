classdef Ctrller < matlab.mixin.Copyable
    % Contains controller implementation matrices (Rc, Mc)
    % Inherits handle class with deep copy functionality
    
    properties
      Rc_;
      Mc_;
      
      solveStatus_; % cvx_status
    end
    
    methods(Static)
      function obj = ctrller_from_cl_maps(clMaps)
        % use closed loop maps as controller
        obj.Rc_ = clMaps.R_;
        obj.Mc_ = clMaps.M_;
        
        obj.solveStatus_ = clMaps.solveStatus_;
      end
    end
    
end