function [RZeros, MZeros] = delay_constraints(sys, Tc, delay)
% TODO: copied from state_fdbk_sls

for t=1:Tc
    RZeros{t} = false(sys.Nx, sys.Nx);
    MZeros{t} = false(sys.Nu, sys.Nx);
    
    for i=1:sys.Nx % note: distance metric only works for chains rn
        for j=1:sys.Nx % for R
            if t < delay * abs(i-j)
                RZeros{t}(i,j) = true;
            end
        end            
    end

    if t > 1 % since Mc{1} == M{1}, don't interfere with that
        MZeros{t} = abs(sys.B2)' * RZeros{t} > 0; 
    end
end
