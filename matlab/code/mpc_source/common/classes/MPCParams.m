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
     
        % Optional params -----------------------------------------
        constrMtx_;
        
        constrUpperbnd_; % upper bound on constrMtx * [state; input]
        constrLowerbnd_; % lower bound on constrMtx * [state; input]
        
        stateUpperbnd_; % upper bound on state (x)
        stateLowerbnd_; % lower bound on state (x)
        
        inputUpperbnd_; % upper bound on input (u)
        inputLowerbnd_; % lower bound on input (u)
        
    end
    
    methods
      function sanity_check_alg_1(obj)
          sanity_check_constr(obj);
          sanity_check_mode(obj);
          
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
              sls_error('One or more required parameters is missing!')
          end
      end
      
      function sanity_check_alg_2(obj)
          sanity_check_alg_1(obj);
          
          e1 = isempty(obj.maxItersCons_);
          e2 = isempty(obj.mu_);
          e3 = isempty(obj.eps_x_);
          e4 = isempty(obj.eps_z_);
          
          if obj.solnMode_ == MPCSolMode.Explicit
              sls_error('Explicit solver mode not available for Algorithm 2!')
          end
          
          if (e1 || e2 || e3 || e4)
              sls_error('One or more required parameters is missing!')
          end 
      end
      
      function sanity_check_constr(obj)
          e1 = isempty(obj.constrMtx_);
          e2 = isempty(obj.constrUpperbnd_);
          e3 = isempty(obj.constrLowerbnd_);
          
          if (~e1 && e2 && e3)
              sls_warning('Constraint matrix specified but no bounds specified! Ignoring')
          elseif (e1 && ~e2 && ~e3)
              sls_warning('Constraint bounds specified but no matrix specified! Ignoring')
          end
      end
      
      function sanity_check_mode(obj)
          e1 = isempty(obj.constrMtx_);
          e2 = isempty(obj.constrUpperbnd_);
          e3 = isempty(obj.constrLowerbnd_);
          e4 = isempty(obj.stateUpperbnd_);
          e5 = isempty(obj.stateLowerbnd_);
          e6 = isempty(obj.inputUpperbnd_);
          e7 = isempty(obj.inputLowerbnd_);
          
          constr = (~e1 && ~e2) || (~e1 && ~e3);
          constr = constr || ~e4 || ~e5 || ~e6 || ~e7;
          
          if obj.solnMode_ == MPCSolMode.UseSolver
              if ~constr
                  sls_warning('No constraints were specified, are you sure you want to use solver?')
                  sls_warning('Closed form would be much, much faster.')
              end              
          elseif obj.solnMode_ == MPCSolMode.ClosedForm
              if constr
                  sls_error('Constraints were specified, cannot use closed form')
              end
          end
      end
    end
end