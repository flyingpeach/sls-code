function [x, u] = simulate_two_blend(sys, ctrller1, ctrller2, simParams, satParams)

simParams.sanity_check();

R1 = ctrller1.Rc_;
M1 = ctrller1.Mc_;
R2 = ctrller2.Rc_;
M2 = ctrller2.Mc_;

T    = length(R1);
tSim = simParams.tSim_;
eta1 = satParams.eta_(1);
eta2 = satParams.eta_(2);

x     = zeros(sys.Nx, tSim); 
u     = zeros(sys.Nu, tSim);
x_hat = zeros(sys.Nx, tSim);
w_hat = zeros(sys.Nx, tSim);

for t=1:1:tSim-1
    for k=1:1:min(t-1, T)
        w_1    = sat(w_hat(:,t-k), eta1);
        w_2    = sat(w_hat(:,t-k), eta2);
        u(:,t) = u(:,t) + M1{k}*w_1 + M2{k}*(w_2 - w_1);
    end

    u(:,t) = sat(u(:,t), satParams.uMax_);
    
    for k=1:1:min(t-1, T-1)
        w_1        = sat(w_hat(:,t-k), eta1);
        w_2        = sat(w_hat(:,t-k), eta2);
        x_hat(:,t) = x_hat(:,t) + R1{k+1}*w_1 + R2{k+1}*(w_2 - w_1);

        if k <= satParams.tau_+1
            x_hat(:,t) = x_hat(:,t) + sys.A^(k-1)*(w_hat(:,t-k) - w_2);
        end
    end

    x(:,t+1)   = sys.A*x(:,t) + sys.B1*simParams.w_(:,t) + sys.B2*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end
