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

        solnMode_; % MPCSolnMode indicating solution mode
        
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
          e8  = isempty(obj.solnMode_);
          e9  = isempty(obj.Q_);
          e10 = isempty(obj.R_);
          
          if (e1 || e2 || e3 || e4 || e5 || e6 || e7 || e8 || e9 || e10)
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
      
      function hasCons = sanity_check_cons(mtx, ub, lb)
          hasMtx = ~isempty(mtx);
          hasUB  = ~isempty(ub);
          hasLB  = ~isempty(lb);

          hasCons = hasMtx;
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
      
      function sanity_check_alg_1(obj)
          if ~isdiag(obj.Q_) || ~isdiag(obj.R_) || ~isdiag(obj.stateConsMtx_) || ~isdiag(obj.inputConsMtx)
              mpc_error('Cannot use Algorithm 1, cost or constraint matrices induce coupling!')
          end
          
          switch obj.solnMode_
              case MPCSolMode.ClosedForm
                  if has_state_cons(obj) || has_input_cons(obj)
                      mpc_error('Cannot use closed form, there are constraints!');
                  end
              case MPCSolMode.Explicit
                  % do nothing
              case MPCSolMode.UseSolver
                  mpc_warning('UseSolver was specified but Explicit would be much faster');
              otherwise
                  mpc_error('Unrecognized solution mode specified');
          end
      end
      
      function sanity_check_alg_2(obj)
          switch obj.solnMode_
              case MPCSolMode.ClosedForm
                  if has_state_cons(obj) || has_input_cons(obj)
                      mpc_error('Cannot use closed form, there are constraints!');
                  end
              case MPCSolMode.Explicit
                  mpc_error('Explicit solutions not available for Algorithm 2!');
              case MPCSolMode.UseSolver
                  % do nothing
              otherwise
                  mpc_error('Unrecognized solution mode specified');
          end         
      end
      
    end
end