classdef CtrllerSettings < matlab.mixin.Copyable
    % Contains settings for alternate implementation finder
    % Note that depending on the mode specified not all params may be used
    % Inherits handle class with deep copy functionality
    
    properties  
      mode_;       % CtrllerMode (i.e. OptL1, OptL1DelayLocal)
      eps_nullsp_; % "nullspace" used in optimization consists of right
                   % singular vectors corresponding to singular values 
                   % below eps_nullsp_

      % used for OptL1DelayLocal mode only
      actDelay_; % actuation delay
      cSpeed_;   % communication speed
      d_;        % d-hop locality constraint
      
      % used for EncourageDelay mode only
      fastCommPen_; % penalize Rc/Mc that would require fast communication
      
      % used for EncourageLocal mode only
      nonLocalPen_; % penalize Rc/Mc with nonlocal behaviour
    end

    methods
      function obj = CtrllerSettings()
        % initialize to zero instead of empty array
        obj.mode_        = 0;
        obj.eps_nullsp_  = 0;
        obj.actDelay_    = 0;
        obj.cSpeed_      = 0;
        obj.d_           = 0;
        obj.fastCommPen_ = 0;
        obj.nonLocalPen_ = 0;
      end      
      
      function statusTxt = sanity_check(obj)
        paramStr = sprintf(', eps_nullsp=%0.2e', obj.eps_nullsp_);
        if not(obj.eps_nullsp_)
            disp('[SLS WARNING] eps_nullsp_=0 causes optimization space to be tiny!');
        end
        
        switch obj.mode_ % check mode & needed params
            case CtrllerMode.OptL1 
                modeStr = 'centralized';            
            case CtrllerMode.OptL1Delayed
                paramStr = check_delay(obj, paramStr);
                modeStr  = 'localized';
            case CtrllerMode.OptL1Localized
                paramStr = check_locality(obj, paramStr);
                modeStr  = 'delayed';
            case CtrllerMode.OptL1DAndL
                paramStr = check_delay(obj, paramStr);
                paramStr = check_locality(obj, paramStr);
                modeStr = 'delayed and localized';
            case CtrllerMode.EncourageDelay
                if not(obj.fastCommPen_)
                    error('[SLS ERROR] Did you forget to specify fastCommPen?');
                end
                modeStr  = 'encourage delay';
                paramStr = [paramStr, sprintf(', fastCommPen=%d', obj.fastCommPen_)];
            case CtrllerMode.EncourageLocal
                if not(obj.nonLocalPen_)
                    error('[SLS ERROR] Did you forget to specify nonLocalPen?');
                end
                modeStr  = 'encourage locality';
                paramStr = [paramStr, sprintf(', nonLocalPen=%d', obj.nonLocalPen_)];
            otherwise
                error('[SLS ERROR] Finding (Rc, Mc) but mode unknown or unspecified!');
        end
        statusTxt = [modeStr, paramStr];
      end
        
      function paramStr = check_locality(obj, paramStr)
        % ensure all needed params for locality exist
        if not(obj.d_)
            disp('[SLS WARNING] Solving with locality constraints but d=0. Did you forget to specify it?');
        end
        paramStr = [paramStr, sprintf(', d=%d', obj.d_)];
      end
      
      function paramStr = check_delay(obj, paramStr)
        % ensure all needed params for delay exist
        if not(obj.cSpeed_)            
            disp('[SLS WARNING] Solving with delay constraints but cSpeed=0. Did you forget to specify it?');
        end
        if not(obj.actDelay_)
            disp('[SLS WARNING] Solving with delay constraints but actDelay=0. Did you forget to specify it?');
        end
        paramStr = [paramStr, sprintf(', cSpeed=%0.2f, actDelay=%d', ...
                                      obj.cSpeed_, obj.actDelay_)];
      end
    end
end