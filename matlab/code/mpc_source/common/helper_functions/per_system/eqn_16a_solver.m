function [Phi_loc, time] = eqn_16a_solver(x_ri, psi_rowi, lamb_rowi, n, mParams)

    rho = mParams.rho_;
    
    model.Q     = sparse(x_ri*x_ri' + rho/2*eye(n));
    model.obj   = rho*(-psi_rowi + lamb_rowi);
    model.A     = sparse([x_ri'; x_ri']);
    model.rhs   = [mParams.state_upperbnd_, mParams.state_lowerbnd_];
    model.sense = '<>';   
    model.lb    = -100*ones(n, 1); % TODO: hacky?
    
    gParams.outputflag = 0;
    result = gurobi(model, gParams);
    time   = result.runtime;
    
    Phi_loc = result.x(:);
end