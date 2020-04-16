classdef MPCParams < matlab.mixin.Copyable
    % Contains parameters for MPC problems
    
    properties
        d_;        % locality
        
        tFIR_;     % order of controller (# spectral components)
        tHorizon_; % MPC time horizon
                
        maxIters_; % maximum allowed iterations for ADMM 
        rho_;      % ADMM update step size
        
        % determines exit conditions at Step 9 (Alg1) / Step 14 (Alg2)
        eps_p_; % convergence criterion for ||Phi(k+1) - Psi(k+1)|| 
        eps_d_; % convergence criterion for ||Psi(k+1) - Psi(k)||
        
        % Algorithm 2 only ----------------------------------------

        maxItersCons_; % maximum allowed iterations for ADMM consensus
        mu_;           % ADMM consensus update step size

        % determines exit conditions at Step 9 (Alg2)
        eps_x_; % convergence criterion for ||X(n+1) - Z(n+1)||
        eps_z_; % convergence criterion for ||Z(n+1) - Z(n)||
    end
end