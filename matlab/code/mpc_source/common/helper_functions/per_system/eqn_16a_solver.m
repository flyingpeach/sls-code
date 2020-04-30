function [Phi_loc] = eqn_16a_solver(x_ri, psi_rowi, lamb_rowi, n, mParams)

    rho = mParams.rho_;
    
    % set up QP
    % minimize Phi*Q*Phi' + obj*Phi'
    % note: constant terms omitted since we are interested in argmin only
    model.Q     = sparse(x_ri*x_ri' + rho/2*eye(n));
    model.obj   = rho*(-psi_rowi + lamb_rowi);
    
    % s.t. x_ri'*Phi < upperbnd
    %      x_ri'*Phi > lowerbnd
    model.A     = sparse([x_ri'; x_ri']);
    model.rhs   = [mParams.stateUpperbnd_, mParams.stateLowerbnd_];
    model.sense = '<>';
    
    % default lower bound is 0; override
    MPC_LB   = -100;
    model.lb = MPC_LB*ones(n, 1);
    
    % solve QP    
    gParams.outputflag = 0;
    result  = gurobi(model, gParams);
    Phi_loc = result.x(:);
end