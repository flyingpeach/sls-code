classdef AltImplSettings < matlab.mixin.Copyable
    % Contains settings for alternate implementation finder
    % Note that depending on the mode specified not all params may be used
    % Inherits handle class with deep copy functionality
    
    properties  
      mode_;      % AltImplMode (i.e. ImplicitOpt, Analytic)
      
      clDiffPen_; % used for ApproxLeaky mode only
                  % penalize deviation from old CL map

      relaxPct_;  % used for ApproxDrop mode only
                  % % of total constraints that are dropped
      
      tol_;       % tolerance used for rank / nullspace calculations
    end
    
    methods
      function obj = AltImplSettings(tol)
        % initialize to zero instead of empty array
        obj.mode_      = 0;
        obj.clDiffPen_ = 0;
        obj.relaxPct_  = 0;
        obj.tol_       = tol;
      end
      
      function statusTxt = sanity_check(obj)
        paramStr = '';
        switch obj.mode_ % check mode & needed params
            case AltImplMode.ImplicitOpt
                modeStr = 'L1-opt/implicit constr';
            case AltImplMode.ExplicitOpt
                modeStr = 'L1-opt/explicit constr';
            case AltImplMode.NullsOpt
                modeStr = 'L1-opt/null space search'
            case AltImplMode.Analytic
                modeStr  = 'analytic';
            case AltImplMode.ApproxDrop
                if not(obj.relaxPct_)
                    error('[SLS ERROR] Finding alt implementation with dropped constraints but relaxPct=0. Did you forget to specify it?');
                end
                modeStr  = 'L1-opt/relaxed constr';
                paramStr = sprintf(', relaxPct=%d', obj.relaxPct_);
            case AltImplMode.ApproxLeaky
                if not(obj.clDiffPen_)
                    error('[SLS ERROR] Finding alt implementation with leaky constraints but clDiffPen=0. Did you forget to specify it?');
                end
                modeStr  = 'L1-opt/leaky constr';
                paramStr = sprintf(', clDiffPen=%d', obj.clDiffPen_);
            otherwise
                error('[SLS ERROR] Finding alt implementation but mode unknown or unspecified!');
        end
               
        statusTxt = [modeStr, paramStr];
      end    
    end
end