function [Rsupport,Msupport,count] = make_d_localized_constraints(A,B,T,d,comms,ta);

Comm_speed =comms; 

Comms_Adj = abs(A)>0;
LocalityR = Comms_Adj^(d-1)>0;
% LocalityR

% for q = 1:N
%     LocalityR(min(2*N,(q-1)*2+1:q*2+2*(d-1)),min(2*N,(q-1)*2+1:q*2)) = [[1,0;1,1];ones(2*(d-1),2)];
% end
if isempty(B)
    B = zeros(size(A,1),1);
end

count = 0;
for t = 1:T
    Rsupport{t} = min(Comms_Adj^(floor(max(0,comms*(t-ta)))),LocalityR)>0;
    Msupport{t} = (abs(B)'*Rsupport{t})>0;
    count = count + sum(sum(Rsupport{t}))+sum(sum(Msupport{t}));
end