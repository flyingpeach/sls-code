classdef MPCParams < matlab.mixin.Copyable
    % Contains parameters for MPC problems
    
    properties
        % Centralized or distributed
        mode_;
        
        % Algorithm 1 and 2 ---------------------------------------       
        locality_;  % locality (d)
        
        tFIR_;      % order of controller (# spectral components)               
        maxIters_;  % maximum allowed iterations for ADMM 
        rho_;       % ADMM update step size
        
        % determines exit conditions at Step 9 (Alg1) / Step 14 (Alg2)
        eps_p_; % convergence criterion for ||Phi(k+1) - Psi(k+1)|| 
        eps_d_; % convergence criterion for ||Psi(k+1) - Psi(k)||
        
        QSqrt_; % penalty for state (square-root)
        RSqrt_; % penalty for input (square-root)

        % Algorithm 2 only ----------------------------------------
        maxItersCons_; % maximum allowed iterations for ADMM consensus
        mu_;           % ADMM consensus update step size

        % determines exit conditions at Step 9 (Alg2)
        eps_x_; % convergence criterion for ||X(n+1) - Z(n+1)||
        eps_z_; % convergence criterion for ||Z(n+1) - Z(n)||
     
        % Optional params for constraints -------------------------
        stateConsMtx_;        
        stateUB_; % stateConsMtx_ * state <= stateUB_
        stateLB_; % stateConsMtx_ * state >= stateLB_
        
        inputConsMtx_;
        inputUB_; % inputConsMtx_ * input <= inputUB_
        inputLB_; % inputConsMtx_ * input >= inputLB_
        
        % Robust MPC parameters -----------------------------------
        distConsMtx_;
        distUB_; % distConsMtx_ * disturbance <= distUB_
        distLB_; % distConsMtx_ * disturbance >= distLB_
        
        % Advanced options ----------------------------------------
        % Leave these blank unless you really know what you're doing
        
        % Which solver to use for main row update
        % No sanity checks accompany this option
        solverMode_; 
        
        % Adaptive ADMM parameter (outer loop only)
        tau_i_;
        tau_d_;
        muAdapt_;
        rhoMax_;
    end

    methods        
      function sanity_check_params_1(obj)
          e1  = isempty(obj.locality_);
          e2  = isempty(obj.tFIR_);
          e3  = isempty(obj.maxIters_);
          e4  = isempty(obj.rho_);
          e5  = isempty(obj.eps_p_);
          e6  = isempty(obj.eps_d_);
          e7  = isempty(obj.QSqrt_);
          e8  = isempty(obj.RSqrt_);
          
          if (e1 || e2 || e3 || e4 || e5 || e6 || e7 || e8)
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
      
      function sanity_check_state_cons(obj)
          hasStateMtx = ~isempty(obj.stateConsMtx_);
          Nx          = size(obj.QSqrt_, 1);

          % This size constraint is enforced for convenient implementation;
          % in theory we can have constraint matrices of any size
          if hasStateMtx && ~isequal(size(obj.stateConsMtx_), [Nx Nx])
              mpc_error('State constraint matrix is wrong size! Expect Nx by Nx');
          end
          
          if hasStateMtx && ~isequal(size(obj.stateUB_), [Nx 1])
              mpc_error('State upper bound is wrong size! Expect Nx by 1');
          end

          if hasStateMtx && ~isequal(size(obj.stateLB_), [Nx 1])
              mpc_error('State upper bound is wrong size! Expect Nx by 1');
          end
      end

      function sanity_check_input_cons(obj)
          hasInputMtx = ~isempty(obj.inputConsMtx_);
          Nu          = size(obj.RSqrt_, 1);
          
          % This size constraint is enforced for convenient implementation;
          % in theory we can have constraint matrices of any size
          if hasInputMtx && ~isequal(size(obj.inputConsMtx_), [Nu Nu])
              mpc_error('Input constraint matrix is wrong size! Expect Nu by Nu');
          end
          
          if hasInputMtx && ~isequal(size(obj.inputUB_), [Nu 1])
              mpc_error('Input upper bound is wrong size! Expect Nu by 1');
          end

          if hasInputMtx && ~isequal(size(obj.inputLB_), [Nu 1])
              mpc_error('Input upper bound is wrong size! Expect Nu by 1');
          end          
      end         
      
      function sanity_check_adaptive(obj)
          hasAdaptive = ~isempty(obj.tau_i_) || ~isempty(obj.tau_d_) || ~isempty(obj.muAdapt_) || ~isempty(obj.rhoMax_);
          emptyParams = isempty(obj.tau_i_) || isempty(obj.tau_d_) || isempty(obj.muAdapt_) || isempty(obj.rhoMax_);
        
          if hasAdaptive
              if emptyParams
                mpc_error('At least one adaptive ADMM parameter was specified but the others were left empty!');
              end
          end
      end
          
      function sanity_check_dist_cons(obj)
          hasDistMtx = ~isempty(obj.distConsMtx_);
          Nx         = size(obj.QSqrt_, 1);
          
          % This size constraint is enforced for convenient implementation;
          % in theory we can have constraint matrices of any size
          if hasDistMtx && ~isequal(size(obj.distConsMtx_), [Nx Nx])
              mpc_error('Disturbance constraint matrix is wrong size! Expect Nx by Nx');
          end
                    
          if hasDistMtx && ~isequal(size(obj.distUB_), [Nx 1])
              mpc_error('Disturbance upper bound is wrong size! Expect Nx by 1');
          end

          if hasDistMtx && ~isequal(size(obj.distLB_), [Nx 1])
              mpc_error('Disturbance lower bound is wrong size! Expect Nx by 1');
          end       
      end
      
      function sanity_check_alg_1(obj)
          sanity_check_params_1(obj);
          sanity_check_state_cons(obj);
          sanity_check_input_cons(obj);
          sanity_check_dist_cons(obj);
          sanity_check_adaptive(obj);
      end
      
      function sanity_check_alg_2(obj)
          sanity_check_params_2(obj);
          sanity_check_state_cons(obj);
          sanity_check_input_cons(obj);
          sanity_check_adaptive(obj);
      end
      
      function sanity_check_cent(obj)
          e1  = isempty(obj.locality_);
          e2  = isempty(obj.tFIR_);
          e3  = isempty(obj.QSqrt_);
          e4  = isempty(obj.RSqrt_);

          if (e1 || e2 || e3 || e4)
              mpc_error('One or more required parameters is missing!')
          end 
          
          sanity_check_state_cons(obj);
          sanity_check_input_cons(obj);
          sanity_check_dist_cons(obj);
      end          
      
      function hasAdaptiveADMM = has_adaptive_admm(obj)
          hasAdaptiveADMM = ~isempty(obj.muAdapt_);
      end
      
      function hasStateCons = has_state_cons(obj)
          hasStateCons = ~isempty(obj.stateConsMtx_);
      end
      
      function hasInputCons = has_input_cons(obj)
          hasInputCons = ~isempty(obj.inputConsMtx_);
      end
      
      function accForDist = accounts_for_disturbance(obj)
          accForDist = ~isempty(obj.distConsMtx_);
      end
      
      function hasCoupling = has_coupling(obj)
          hasObjCoupling  = ~(isdiag(obj.QSqrt_) && isdiag(obj.RSqrt_));
          hasConsCoupling = ~(isdiag(obj.stateConsMtx_) && isdiag(obj.inputConsMtx_));
          hasDistCoupling = ~isdiag(obj.distConsMtx_);
          
          hasCoupling = hasObjCoupling || hasConsCoupling || hasDistCoupling;          
      end
      
    end
end