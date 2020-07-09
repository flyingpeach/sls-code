classdef SLSParams < matlab.mixin.Copyable
    % Contains parameters for SLS 
    % Inherits handle class with deep copy functionality

    properties
      T_;           % finite impulse response horizon
      constraints_; % list of constraints
      objectives_;  % list of objectives
      approx_;      % boolean indicating whether to relax CL constraints
      stabCoeff_;   % regularization coefficient for stability objective
    end

    methods
      function obj = SLSParams()
        obj.T_           = 0;
        obj.constraints_ = {};
        obj.objectives_  = {};
        obj.approx_      = false;
        obj.stabCoeff_ = 0;
      end

      function add_constraint(obj, type, val)
          % val : value (i.e. locality) associated with constraint
          if not(is_instance(type, 'SLSConstraint'))
              sls_error('Tried to add invalid constraint')
          end
          obj.constraints_ = obj.add(obj.constraints_, type, val);
      end
      
      function add_objective(obj, type, reg_val)
          % reg_val : regularization value for objective
          if not(is_instance(type, 'SLSObjective'))
              sls_error('Tried to add invalid objective')
          end
          obj.objectives_ = obj.add(obj.objectives_, type, reg_val);
      end
            
      function print(obj)
          fprintf(['T = ', num2str(obj.T_), '\n']);
          fprintf('Objectives:\n');  obj.print_(obj.objectives_, ', reg coeff:  ');
          fprintf('Constraints:\n'); obj.print_(obj.constraints_, ': ');
          if obj.approx_
              fprintf(['Using approximate SLS algorithm, coeff = ', ...
                        num2str(obj.stabCoeff_), '\n'])
          end
          fprintf('\n')
      end
      
      function sanity_check(obj)
          if isempty(obj.objectives_)
              sls_warning('No objectives were specified. Only finding a feasible solution');
          end
          
          if obj.approx_ && ~obj.stabCoeff_
              sls_warning('Using approximate SLS but stabCoeff=0; no regularization for stability');
              sls_warning('This may result in an internally unstable controller!');
          end
          
          if obj.stabCoeff_ && ~obj.approx_
              sls_warning('stabCoeff specified but not used');
          end
      end
    end
    
    methods (Static)
      function list = add(list, type, val)
          add_idx = length(list) + 1;
          % if type already exists in list, replace value
          for i=1:length(list)
              thisItem = list{i};
              thisType = thisItem{1};
              if thisType == type
                  sls_warning('The constraint type you tried to add already existed; replaced old value with new value');
                  add_idx = i;
                  break;
              end
          end
          list{add_idx} = {type, val};
      end
      
      function print_(list, linkStr)
          for i=1:length(list)
              thisItem = list{i};
              thisType = thisItem{1};
              thisVal  = thisItem{2};
              fprintf(['\t', char(thisType), linkStr, num2str(thisVal), '\n']);
          end
      end
    end
end