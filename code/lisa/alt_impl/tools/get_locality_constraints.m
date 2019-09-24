function [RSupp, MSupp, count] = get_locality_constraints(sys, Tc, locality)
% Gets supports on R, M based on
% locality. Locality = how many neighbours we communicate with (we will
% NEVER communicate to neighbours outside this area

commsAdj  = abs(sys.A) > 0;
RSupp    = commsAdj^(locality) > 0;
MSupp    = (abs(sys.B2)' * RSupp) > 0; 

count = Tc * (sum(sum(RSupp)) + sum(sum(MSupp)));
end
