function [Phi_loc] = eqn_16a_explicit(x_ri, psi_rowi, lamb_rowi, n, mParams)

    rho = mParams.rho_;
    M   = inv(2*x_ri*x_ri' + rho*eye(n));
    a   = psi_rowi - lamb_rowi;
    
    b1 = mParams.stateUpperbnd_; 
    b2 = mParams.stateLowerbnd_;

    crit  = rho*a*M*x_ri;

    lamb = 0;
    if crit - b1 > 0
        lamb = (crit - b1) / (x_ri'*M*x_ri);
    elseif crit - b2 < 0
        lamb = (crit - b2) / (x_ri'*M*x_ri);
    end

    Phi_loc = (rho*a - lamb*x_ri')*M;
end