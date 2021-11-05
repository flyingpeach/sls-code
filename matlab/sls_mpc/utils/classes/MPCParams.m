classdef MPCParams < matlab.mixin.Copyable
    % Contains parameters for MPC problems
    
    properties
        % Centralized or distributed
        mode_;
                
        % Main parameters -----------------------------------------       
        locality_;  % locality (d)
        
        tFIR_;      % order of controller (# spectral components)               
        maxIters_;  % maximum allowed iterations for ADMM 
        rho_;       % ADMM update step size (initial rho if using adaptive)
        
        % Adaptive ADMM parameter (outer loop only)
        % If primal residue > muAdapt * dual residue,   multiply rho by tau_i
        % If   dual residue > muAdapt * primal residue  divide rho by tau_d
        tau_i_;
        tau_d_;    
        muAdapt_;  
        rhoMax_;   % Maximum value of rho that can be used
        
        % determines exit conditions at Step 9 (Alg1) / Step 14 (Alg2)
        eps_p_; % convergence criterion for ||Phi(k+1) - Psi(k+1)|| 
        eps_d_; % convergence criterion for ||Psi(k+1) - Psi(k)||
        
        QSqrt_; % penalty for state (square-root)
        RSqrt_; % penalty for input (square-root)

        % Consensus parameters ------------------------------------
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
        
        locNoiseBound_; % ||local disturbance||_2 <= locNoiseBound

        % Terminal set parameters ---------------------------------
        % terminal set is specified by terminal_H_ * x <= terminal_h_
        terminal_H_;
        terminal_h_;
        terminal_cost_ = false; % Whether to add terminal cost as well
        
        % Advanced options ----------------------------------------
        % Leave these blank unless you really know what you're doing
        
        % Whether to use solver (instead of explicit) for row update
        useSolver_ = false;        
    end
    
    methods(Static)
      function check_constraint(mtx, ub, lb, msize)
          if ~isequal(size(mtx), [msize msize])
              mpc_error('State/Input/Dist matrix is wrong size!');
          end
          if ~isequal(size(ub), [msize 1])
              mpc_error('State/Input/Dist upper bound is wrong size!');
          end
          if ~isequal(size(lb), [msize 1])
              mpc_error('State/Input/Dist lower bound is wrong size!');
          end
      end
    end
    
    methods
        
      %% Main methods to be called externally
      function sanity_check_dist(obj)
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
          
          check_constraints(obj);
          check_adaptive(obj);
          check_consensus_params(obj);   
          check_terminal(obj);
      end
            
      function sanity_check_cent(obj)
          e1  = isempty(obj.locality_);
          e2  = isempty(obj.tFIR_);
          e3  = isempty(obj.QSqrt_);
          e4  = isempty(obj.RSqrt_);

          if (e1 || e2 || e3 || e4)
              mpc_error('One or more required parameters is missing!')
          end 
          
          check_constraints(obj);
          check_terminal(obj);          
      end       
      
      %% Helpers      
      function check_consensus_params(obj)
          if need_consensus(obj)
              e1 = isempty(obj.maxItersCons_);
              e2 = isempty(obj.mu_);
              e3 = isempty(obj.eps_x_);
              e4 = isempty(obj.eps_z_);          

              if (e1 || e2 || e3 || e4)
                  mpc_error('One or more required consensus parameters is missing!')
              end
          end
      end
      
      function check_adaptive(obj)
          if has_adaptive_admm(obj)          
              emptyParams = isempty(obj.tau_i_) || isempty(obj.tau_d_) || isempty(obj.muAdapt_) || isempty(obj.rhoMax_);
              if emptyParams
                mpc_error('At least one adaptive ADMM parameter was specified but the others were left empty!');
              end
          end
      end
      
      function check_constraints(obj)
          Nx = size(obj.QSqrt_, 1);
          Nu = size(obj.RSqrt_, 1);
          if has_state_cons(obj)
              MPCParams.check_constraint(obj.stateConsMtx_, obj.stateUB_, obj.stateLB_, Nx);
          end
          if has_input_cons(obj)
              MPCParams.check_constraint(obj.inputConsMtx_, obj.inputUB_, obj.inputLB_, Nu);
          end
          if has_polytopic_noise(obj)
              MPCParams.check_constraint(obj.distConsMtx_, obj.distUB_, obj.distLB_, Nx);
          end
      end
      
      function check_terminal(obj)
          if has_terminal_set(obj) || obj.terminal_cost_
              emptyParams = isempty(obj.terminal_H_) || isempty(obj.terminal_h_);
              if emptyParams
                mpc_error('At least one terminal parameter was specified but the others were left empty!');
              end              
          end          
      end
            
      %% Boolean functions
      function hasStateCons = has_state_cons(obj)
          hasStateCons = ~isempty(obj.stateConsMtx_) || ~isempty(obj.stateUB_) || ~isempty(obj.stateLB_);
      end
      
      function hasInputCons = has_input_cons(obj)
          hasInputCons = ~isempty(obj.inputConsMtx_) || ~isempty(obj.inputUB_) || ~isempty(obj.inputLB_);
      end
      
      function hasCoupling = has_coupling(obj)
          hasObjCoupling  = ~(isdiag(obj.QSqrt_) && isdiag(obj.RSqrt_));
          hasConsCoupling = ~(isdiag(obj.stateConsMtx_) && isdiag(obj.inputConsMtx_));
          hasDistCoupling = ~isdiag(obj.distConsMtx_);
          
          hasCoupling = hasObjCoupling || hasConsCoupling || hasDistCoupling;          
      end

      function hasTerminalSet = has_terminal_set(obj)
          hasTerminalSet = ~isempty(obj.terminal_H_) || ~isempty(obj.terminal_h_);
      end
      
      function hasAdaptiveADMM = has_adaptive_admm(obj)
          hasAdaptiveADMM = ~isempty(obj.tau_i_) || ~isempty(obj.tau_d_) || ~isempty(obj.muAdapt_) || ~isempty(obj.rhoMax_);
      end

      function needConsensus = need_consensus(obj)
          % Only the coupled, non-robust case needs consensus
          needConsensus = has_coupling(obj) && ~is_robust(obj);
      end

      function isRobust = is_robust(obj)
          isRobust = has_polytopic_noise(obj) || has_loc_bounded_noise(obj);
      end

      function hasPolyNoise = has_polytopic_noise(obj)
          hasPolyNoise = ~isempty(obj.distConsMtx_) || ~isempty(obj.distUB_) || ~isempty(obj.distLB_);
      end

      function hasLocBoundNoise = has_loc_bounded_noise(obj)
          hasLocBoundNoise = ~isempty(obj.locNoiseBound_);
      end      
      
      
    end
end