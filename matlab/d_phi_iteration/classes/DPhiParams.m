classdef DPhiParams < matlab.mixin.Copyable
    % Contains parameters for simulating the system
    % Inherits handle class with deep copy functionality

    properties
      % Core SLS properties
      locality_; % size of locality, 'd' in papers; 1 means self-comm only
      tFIR_;     % FIR length of Phi
      
      % Assume that regulated output z takes the form z = Hx*x + Hu*u
      Hx_;
      Hu_;
      
      % Assume that performance objective is x'*Qx + u'*Ru
      % Where Q = QSqrt'*QSqrt and R = RSqrt'*RSqrt
      QSqrt_;
      RSqrt_;
      
      % Parameters for algorithm       
      stabType_; % StabType.L1, .LInf, or .Nu
      betaStep_   = 5e-2;  % Minimum robust bound progress per step
      betaStop_   = -inf;  % Stop when we reach this beta or when infeasible
      randomizeD_ = false; % Default is to minimize D in D-step      
    end
    
    methods
      function sanity_check(obj)
          e1  = isempty(obj.locality_);
          e2  = isempty(obj.tFIR_);
          e3  = isempty(obj.Hx_);
          e4  = isempty(obj.Hu_);
          e5  = isempty(obj.QSqrt_);
          e6  = isempty(obj.RSqrt_);
          e7  = isempty(obj.stabType_);
          
          if (e1 || e2 || e3 || e4 || e5 || e6 || e7)
              error('One or more required parameters is missing!\n')
          end          
      end
    end
    
end