classdef PSeq_Params  < handle
    % Hold paramaters that will be global for all sequence elements.

    properties
        sys;
        channels = {'x', 'y', 'z'};
        rf_spoil = true;
        rf_spoil_idx = 0;
        rf_spoil_phase = 0;

    end

    methods
        function obj = PSeq_Params(varargin)
            % Construct an instance of this class
            %   
            % Params
            % ---------
            % specs is a class with MaxGrad and MaxSlew 
            %    TODO: Switch to inputparser for this

            p = inputParser;
            addParameter(p, 'max_grad', 30);
            addParameter(p, 'max_slew', 100);
            addParameter(p, 'grad_raster', 10e-6);

            parse(p, varargin{:});
            inputs = p.Results;

            obj.sys = mr.opts('MaxGrad',inputs.max_grad, ...
                              'GradUnit','mT/m',...
                              'MaxSlew',inputs.max_slew, ...
                              'SlewUnit','T/m/s',...
                              'rfRingdownTime', 60e-6, ...
                              'rfDeadtime', 100e-6, ...
                              'adcDeadTime', 40e-6, ...
                              'B0', 2.98, ...
                              'adcSamplesLimit',8192, ...
                              'adcRasterTime', 2.0e-6, ...
                              'gradRasterTime', inputs.grad_raster, ...
                              'blockDurationRaster', inputs.grad_raster);  
        end

        function increment_rf_spoiling(obj)
            % Increment the rf spoiling phase.
            % 
            % rf_spoil_phase should then be used any time an RF or ADC even is played.
            % 
            % TODO
            % ----
            % * Add support for different style of rf spoiling (mainly random, or list based)

            obj.rf_spoil_idx = obj.rf_spoil_idx + 1;
            obj.rf_spoil_phase = (117 * pi / 180) * obj.rf_spoil_idx * obj.rf_spoil_idx / 2;
            obj.rf_spoil_phase = mod(obj.rf_spoil_phase, 2*pi);
        end
    end
end