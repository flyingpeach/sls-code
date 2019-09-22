function [RZeros, MZeros] = get_delay_constraints(sys, Tc, delay)
% Gets zero-constraints on R, M (i.e. where they must be zero) based on
% delay specified. Delay = amount of time it takes to go from one node
% to its neighbour

for t=1:Tc
    RZeros{t} = false(sys.Nx, sys.Nx);    
    for i=1:sys.Nx % note: distance metric only works for chains rn
        for j=1:sys.Nx
            if t < delay * abs(i-j)
                RZeros{t}(i,j) = true;
            end
        end            
    end
    MZeros{t} = abs(sys.B2)' * RZeros{t} > 0; 
end
