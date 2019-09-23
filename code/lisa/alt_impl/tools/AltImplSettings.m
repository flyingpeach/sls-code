classdef AltImplSettings < matlab.mixin.Copyable
    % Contains settings for alternate implementation finder
    % Note that depending on the mode specified not all params may be used
    % Inherits handle class with deep copy functionality
    
    properties  
      mode_;      % AltImplMode (i.e. ExactOpt, Analytic)

      tol_;       % tolerance used in non-relaxed rank/nullsp calculations
      
      svThresh_;  % used in ApproxLS mode
                  % use right singular vectors corresponding to singular 
                  % values below this threshold

      delay_;     % used in StrictDelay mode
                  % amount of delay to enforce for strict locality

      m1NonzeroPen_ % used in StrictDelay mode
                    % Mc1 = M1 is requirement; penalize nonzero values that
                    % violate delay constraints
                  
      clDiffPen_; % used in ApproxLeaky mode
                  % penalize deviation from old CL map
      
      fastCommPen_; % used in EncourageDelay mode
                    % penalize Rc/Mc that would require fast communication
    end
    
    methods
      function obj = AltImplSettings()
        % initialize to zero instead of empty array
        obj.mode_        = 0;
        
        obj.tol_         = 0;
        obj.svThresh_    = 0;
        obj.delay_       = 0;
        obj.clDiffPen_   = 0;
        obj.fastCommPen_ = 0;
      end      
      
      function statusTxt = sanity_check(obj)
        paramStr = '';
     
        switch obj.mode_ % check mode & needed params
            case AltImplMode.ExactOpt
                if not(obj.tol_)
                    disp('[SLS WARNING] Zero tolerance may make feasible problems seem infeasible!');
                end                
                modeStr = 'exact constraints';
                paramStr = sprintf(', tol=%0.2e', obj.tol_);
            case AltImplMode.Analytic
                modeStr  = 'analytic';
            case AltImplMode.ApproxLS
                if not(obj.svThresh_)
                    error('[SLS ERROR] Did you forget to specify svThresh?');
                end
                modeStr = 'relaxed nullspace';
                paramStr = sprintf(', svThresh=%0.2e', obj.svThresh_);
            case AltImplMode.ApproxLeaky
                if not(obj.clDiffPen_)
                    error('[SLS ERROR] Did you forget to specify clDiffPen?');
                end
                if not(obj.tol_)
                    disp('[SLS WARNING] Zero tolerance may make feasible problems seem infeasible!');
                end                
                modeStr  = 'leaky constr';
                paramStr = sprintf(', tol=%0.2e, clDiffPen=%d', obj.clDiffPen_, obj.tol_);
            case AltImplMode.StrictDelay
                if not(obj.tol_)
                    disp('[SLS WARNING] Zero tolerance may make feasible problems seem infeasible!');
                end                
                if not(obj.m1NonzeroPen_)
                    error('[SLS ERROR] Did you forget to specify m1NonzeroPen?');
                end
                if not(obj.delay_)
                    error('[SLS ERROR] Did you forget to specify delay?');
                end
                modeStr = 'strict delay';
                paramStr = sprintf(', m1NonzeroPen=%0.2e, delay=%0.1f', obj.m1NonzeroPen_, obj.delay_);
            case AltImplMode.EncourageDelay
                if not(obj.tol_)
                    disp('[SLS WARNING] Zero tolerance may make feasible problems seem infeasible!');
                end                
                if not(obj.fastCommPen_)
                    error('[SLS ERROR] Did you forget to specify fastCommPen?');
                end
                modeStr  = 'encourage delay';
                paramStr = sprintf(', tol=%0.2e, fastCommPen=%d', obj.clDiffPen_, obj.fastCommPen_);
            
            otherwise
                error('[SLS ERROR] Finding alt implementation but mode unknown or unspecified!');
        end
               
        statusTxt = [modeStr, paramStr];
      end    
    end
end