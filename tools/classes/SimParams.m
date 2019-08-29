classdef SimParams < matlab.mixin.Copyable
    % inherits handle class with deep copy functionality
    % contains parameters for simulating the system

    properties  
      tSim_;     % amount of time to simulate
      w_;        % disturbance w(t) (as per (3.1)
      openLoop_; % whether to simulate the system open loop only
    end
end