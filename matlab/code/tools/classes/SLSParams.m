classdef SLSParams < matlab.mixin.Copyable
    % Contains parameters for SLS 
    % Inherits handle class with deep copy functionality

    properties
      T_;           % finite impulse response horizon
      constraints_; % list of constraints
      objectives_;  % list of objectives
      approx_;      % boolean indicating whether to relax CL constraints
      approxCoeff_; % regularization coefficient for stability objective
    end

    methods
      function obj = SLSParams()
        obj.T_           = 0;
        obj.constraints_ = {};
        obj.objectives_  = {};
        obj.approx_      = false;
        obj.approxCoeff_ = 0;
      end
      
      function obj = add_constraint(obj, type, val)
          % val : value (i.e. locality) associated with constraint
          %       unused if type = Standard or Robust
          if not(is_instance(type, 'SLSConstraint'))
              sls_error('Tried to add invalid constraint')
          end
          obj.constraints_{length(obj.constraints_)+1} = [type, val];
      end
      
      function obj = add_objective(obj, type, reg_val)
          % reg_val : regularization value for objective
          if not(is_instance(type, 'SLSObjective'))
              sls_error('Tried to add invalid objective')
          end
          obj.objectives_{length(obj.objectives_)+1} = [type, reg_val];
          
      end
      
      function print_constraints(obj)
          for i=1:length(obj.constraints_)
              this_constr = obj.constraints_{i};
              constrType  = this_constr(1);
              constrVal   = this_constr(2);
              fprintf([class(constrType), ': ', num2str(constrVal), '\n']);
          end
      end
      
      function print_objectives(obj)
          for i=1:length(obj.objectives_)
              this_obj  = obj.objectives_{i};
              objType   = this_obj(1);
              objRegVal = this_obj(2);
              fprintf([class(objType), ', reg const: ' , num2str(objRegVal), '\n']);
          end
      end
      
      function print(obj)
          fprintf(['T = ,' obj.T_, '\n']);
          fprintf('Objectives:\n');  print_objectives(obj);
          fprintf('Constraints:\n'); print_constraints(obj);
          if obj.approx_
              fprintf('Using approximate SLS algorithm')
          end          
      end
      
      function sanity_check(obj)
          if isempty(obj.objectives_)
              sls_warning('No objectives were specified. Only finding a feasible solution');
          end
          
          if obj.approx_ && ~obj.approxCoeff_
              sls_warning('Using approximate SLS but approxCoeff=0; no regularization for stability');
              sls_warning('This may result in an internally unstable controller!');
          end
          
          if obj.approxCoeff_ && ~obj.approx_
              sls_warning('approxCoeff specified but not used');
          end
      end
    end
end