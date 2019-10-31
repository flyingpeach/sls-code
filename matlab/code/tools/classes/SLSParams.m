classdef SLSParams < matlab.mixin.Copyable
    % Contains parameters for SLS 
    % Note that depending on the solver called, not all params may be used
    % Inherits handle class with deep copy functionality
    
    properties  
      mode_;     % an SLSMode (i.e. Basic, DLocalized)
      rfd_;      % whether we want rfd
      
      % used for all modes
      tFIR_;     % finite impulse response horizon
      
      % used for DLocalized and ApproxDLocalized modes only
      actDelay_; % actuation delay
      cSpeed_;   % communication speed
      d_;        % d-hop locality constraint
      
      robCoeff_; % regularization coeff for robust stability (used in approx-d sls only)
      rfdCoeff_; % regularization coeff for rfd (used in rfd only)
      
      obj_;      % an Objective (i.e. H2, HInf)
    end
    
    methods
      function obj = SLSParams()
        % initialize to zero instead of empty array
        obj.mode_     = 0; 
        obj.rfd_      = 0;
        
        obj.tFIR_     = 0;        
        obj.actDelay_ = 0;
        obj.cSpeed_   = 0;
        obj.d_        = 0;
        
        obj.robCoeff_ = 0; 
        obj.rfdCoeff_ = 0;
        
        obj.obj_      = 0;
      end
      
      function statusTxt = sanity_check(obj)
        if not(obj.tFIR_)
            error('[SLS ERROR] tFIR=0. Did you forget to specify it?');
        end

        paramStr = [char(10), char(9), sprintf('tFIR=%d', obj.tFIR_)];
        
        switch obj.mode_ % check mode & needed params
            case SLSMode.Basic
                modeStr = 'basic';
            case SLSMode.DLocalized
                paramStr = check_d_local(obj, paramStr);
                modeStr = 'd-localized';
            case SLSMode.ApproxDLocalized
                if not(obj.robCoeff_)
                    error('[SLS ERROR] Solving with approximate d-localized SLS but robCoeff=0. Did you forget to specify it?');
                end
                paramStr = check_d_local(obj, paramStr);
                paramStr = [paramStr, sprintf(', robCoeff=%0.2f', obj.robCoeff_)];
                modeStr  = 'approx d-localized';
            otherwise
                error('[SLS ERROR] SLS mode unknown or unspecified!');
        end

        switch obj.obj_ % check objective & needed params
            case Objective.H2
                objStr = 'H2';
            case Objective.HInf
                objStr = 'HInf';
            case Objective.L1
                objStr = 'L1';
            otherwise
                objStr = 'constant';
        end
        
        statusTxt = [modeStr, ' SLS with ', objStr, ' objective'];

        if obj.rfd_ == true
            if not(obj.rfdCoeff_)
                disp('[SLS WARNING] Solving with RFD but rfdCoeff=0. Did you forget to specify it?');
            end
            statusTxt = [statusTxt, ' and RFD, ', sprintf('rfdCoeff=%0.2f', obj.rfdCoeff_)];
        end
        statusTxt = [statusTxt, paramStr];
      end

      function paramStr = check_d_local(obj, paramStr)
        % ensure all the needed parameters for localized sls are specified
        if not(obj.d_)
            disp('[SLS WARNING] Solving with locality constraints but d=0. Did you forget to specify it?');
        end
        if not(obj.cSpeed_)            
            disp('[SLS WARNING] Solving with locality constraints but cSpeed=0. Did you forget to specify it?');
        end
        if not(obj.actDelay_)
            disp('[SLS WARNING] Solving with locality constraints but actDelay=0. Did you forget to specify it?');
        end
        
        paramStr = [paramStr, sprintf(', d=%d, cSpeed=%0.2f, actDelay=%d', ...
                              obj.d_, obj.cSpeed_, obj.actDelay_)];
      end
    end
    
end