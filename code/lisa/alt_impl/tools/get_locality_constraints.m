function [RZeros, MZeros] = get_locality_constraints(sys, locality)
% Gets zero-constraints on R, M (i.e. where they must be zero) based on
% locality. Locality = how many neighbours we communicate with (we will
% NEVER communicate to neighbours outside this area

commsAdj  = abs(sys.A) > 0;
RZeros    = not(commsAdj^(locality) > 0);
MZeros    = (abs(sys.B2)' * RZeros) > 0; 
end
