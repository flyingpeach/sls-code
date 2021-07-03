classdef MPCWarmStart < matlab.mixin.Copyable
    % Matrices to be passed on to the next timestep for a warm start
    
    properties
        Psi_;
        Lambda_;
    end

end