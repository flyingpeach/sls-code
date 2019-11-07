classdef CtrllerParams < matlab.mixin.Copyable
    % Contains parameters for controller implementation matrices
    % Note that depending on the mode specified not all params may be used
    % Inherits handle class with deep copy functionality
    
    properties  
      mode_;       % SLSMode (i.e. Basic, Localized)
      
      % used for all modes
      eps_nullsp_; % "nullspace" used in optimization consists of right
                   % singular vectors corresponding to singular values 
                   % below eps_nullsp_

      tc_;         % number of spectral components for Rc, Mc
                   
      % used for Delayed / Localized / DAndL / ApproxDAndL modes
      actDelay_; % actuation delay
      cSpeed_;   % communication speed
      d_;        % d-hop locality constraint
      
      % used for EncourageDelay mode only
      fastCommPen_; % penalize Rc/Mc that would require fast communication
      
      % used for EncourageLocal mode only
      nonLocalPen_; % penalize Rc/Mc with nonlocal behaviour
    end

    methods
      function obj = CtrllerParams()
        % initialize to zero instead of empty array
        obj.mode_        = 0;
        obj.eps_nullsp_  = 0;
        obj.tc_          = 0;
        obj.actDelay_    = 0;
        obj.cSpeed_      = 0;
        obj.d_           = 0;
        obj.fastCommPen_ = 0;
        obj.nonLocalPen_ = 0;
      end      
      
      function statusTxt = sanity_check(obj)
        if not(obj.tc_)
            error('[SLS ERROR] tc=0. Did you forget to specify it?');
        end
        if not(obj.eps_nullsp_)
            disp('[SLS WARNING] eps_nullsp=0 causes optimization space to be tiny!');
        end
        paramStr = [char(10), char(9), sprintf(', tc=%d, eps_nullsp=%0.2e', obj.tc_, obj.eps_nullsp_)];
        
        switch obj.mode_ % check mode & needed params
            case SLSMode.Basic 
                modeStr = 'centralized';            
            case SLSMode.Delayed
                paramStr = check_delay(obj, paramStr);
                modeStr  = 'localized';
            case SLSMode.Localized
                paramStr = check_locality(obj, paramStr);
                modeStr  = 'delayed';
            case SLSMode.DAndL
                paramStr = check_delay(obj, paramStr);
                paramStr = check_locality(obj, paramStr);
                modeStr = 'delayed and localized';
            case SLSMode.EncourageDelay
                if not(obj.fastCommPen_)
                    error('[SLS ERROR] Did you forget to specify fastCommPen?');
                end
                modeStr  = 'encourage delay';
                paramStr = [paramStr, sprintf(', fastCommPen=%d', obj.fastCommPen_)];
            case SLSMode.EncourageLocal
                if not(obj.nonLocalPen_)
                    error('[SLS ERROR] Did you forget to specify nonLocalPen?');
                end
                modeStr  = 'encourage locality';
                paramStr = [paramStr, sprintf(', nonLocalPen=%d', obj.nonLocalPen_)];
            otherwise
                errStr = ['An invalid SLSMode was chosen for CtrllerParams!', char(10), ...
                          'Options: Basic, Delayed, Localized, DAndL, EncourageDelay, EncourageLocal'];
                error(errStr);
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