classdef MPCParams < matlab.mixin.Copyable
    % Contains parameters for MPC problems
    
    properties        
        % Algorithm 1 and 2 ---------------------------------------       
        locality_;  % locality (d)
        
        tFIR_;      % order of controller (# spectral components)
        tHorizon_;  % MPC time horizon
                
        maxIters_;  % maximum allowed iterations for ADMM 
        rho_;       % ADMM update step size
        
        % determines exit conditions at Step 9 (Alg1) / Step 14 (Alg2)
        eps_p_; % convergence criterion for ||Phi(k+1) - Psi(k+1)|| 
        eps_d_; % convergence criterion for ||Psi(k+1) - Psi(k)||
        
        Q_; % penalty for state
        R_; % penalty for input

        % Algorithm 2 only ----------------------------------------
        maxItersCons_; % maximum allowed iterations for ADMM consensus
        mu_;           % ADMM consensus update step size

        % determines exit conditions at Step 9 (Alg2)
        eps_x_; % convergence criterion for ||X(n+1) - Z(n+1)||
        eps_z_; % convergence criterion for ||Z(n+1) - Z(n)||
     
        % Optional params for constraints -------------------------
        stateConsMtx_;
        inputConsMtx_;
        
        stateUB_; % upper bound on stateConsMtx_ * state;
        stateLB_; % lower bound on stateConsMtx_ * state;
        
        inputUB_; % upper bound on inputConsMtx_ * input;
        inputLB_; % lower bound on inputConsMtx_ * input;
    end
    
    methods
      function sanity_check_params_1(obj)
          e1  = isempty(obj.locality_);
          e2  = isempty(obj.tFIR_);
          e3  = isempty(obj.tHorizon_);
          e4  = isempty(obj.maxIters_);
          e5  = isempty(obj.rho_);
          e6  = isempty(obj.eps_p_);
          e7  = isempty(obj.eps_d_);
          e8  = isempty(obj.Q_);
          e9  = isempty(obj.R_);
          
          if (e1 || e2 || e3 || e4 || e5 || e6 || e7 || e8 || e9)
              mpc_error('One or more required parameters is missing!')
          end
      end
      
      function sanity_check_params_2(obj)
          sanity_check_params_1(obj);
          
          e1 = isempty(obj.maxItersCons_);
          e2 = isempty(obj.mu_);
          e3 = isempty(obj.eps_x_);
          e4 = isempty(obj.eps_z_);          
          
          if (e1 || e2 || e3 || e4)
              mpc_error('One or more required parameters is missing!')
          end 

      end
      
      function hasthisCons = sanity_check_cons(mtx, ub, lb)
          hasMtx = ~isempty(mtx);
          hasUB  = ~isempty(ub);
          hasLB  = ~isempty(lb);

          hasthisCons = hasMtx;
          if hasMtx && ~hasUB && ~hasLB
              mpc_error('A constraint matrix was specified with no corresponding bounds!');
          elseif ~hasMtx && (hasUB || hasLB)
              mpc_error('Constraint bounds were specified with no corresponding matrix!');
          end
      end
      
      function hasStateCons = has_state_cons(obj)
          hasStateCons = sanity_check_cons(obj.stateConsMtx_, obj.stateUB_, obj.stateLB_);
      end

      function hasInputCons = has_input_cons(obj)
          hasInputCons = sanity_check_cons(obj.inputConsMtx_, obj.inputUB_, obj.inputLB_);
      end
      
      function hasCoupling = has_coupling(obj)
          % at least one non-diagonal cost / constraint matrix
          hasCoupling = ~(isdiag(obj.Q_) && isdiag(obj.R_) && isdiag(obj.stateConsMtx_) && isdiag(obj.inputConsMtx));
      end
      
    end
end