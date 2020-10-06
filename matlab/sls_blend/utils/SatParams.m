classdef SatParams < matlab.mixin.Copyable
    % Contains parameters for SLS with saturation as per 
    % https://arxiv.org/pdf/2006.12766.pdf    
    % Inherits handle class with deep copy functionality

    properties
        noiseStd_; % expected noise
        eta_;      % blending parameter
        tau_;      % anti-windup parameter, used in simulation only
        xMax_;     % maximum allowed state magnitude
        uMax_;     % maximum allowed input magnitude
    end
end