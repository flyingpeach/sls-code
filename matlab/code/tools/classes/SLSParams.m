classdef SLSParams < matlab.mixin.Copyable
    % Contains parameters for SLS 
    % Inherits handle class with deep copy functionality

    properties
      T_;           % finite impulse response horizon
      constraints_; % list of constraints
      objectives_;  % list of objectives
    end

    methods
      function obj = SLSParams()
        obj.T_           = 0;
        obj.constraints_ = {};
        obj.objectives_  = {};
      end
   
      function obj = add_constraint(obj, type, val)
          % val : value (i.e. locality) associated with constraint
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
      
      function print_params(obj)
          fprintf(['T = ,' obj.T_, '\n']);
          fprintf('Objectives:\n');  print_objectives(obj);
          fprintf('Constraints:\n'); print_constraints(obj);
      end
      
    end
end