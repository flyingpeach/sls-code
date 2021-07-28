classdef MPCStats < matlab.mixin.Copyable
    % Contains statistics for DMPC problems: runtime, iters, etc.

    properties
        time_;
        iters_;
        consIters_;
    end

end