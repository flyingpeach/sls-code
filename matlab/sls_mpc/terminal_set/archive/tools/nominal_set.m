function [H_, h_] = nominal_set(Hxu, ZAB, hxu, Nx)
% Gives N-step controllable set for some set defined by
% Hxu*Phi*x <= hxu, ZAB*Phi = [I;0]
% Output set is defined as x: H_*x <= h_

ZABp = pinv(ZAB);
nPhi = size(ZAB, 2);
Q    = Hxu * (eye(nPhi) - ZABp * ZAB);

[H_, h_] = project_set([Q Hxu*ZABp(:,1:Nx)], hxu, Nx);

end