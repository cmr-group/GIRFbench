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
        function obj = PSeq_Params(specs)
            % Construct an instance of this class
            %   
            % Params
            % ---------
            % specs is a class with MaxGrad and MaxSlew 
            %    TODO: Switch to inputparser for this

            obj.sys = mr.opts('MaxGrad',specs.MaxGrad, ...
                              'GradUnit','mT/m',...
                              'MaxSlew',specs.MaxSlew, ...
                              'SlewUnit','T/m/s',...
                              'rfRingdownTime', 30e-6, ...
                              'rfDeadtime', 100e-6, ...
                              'adcDeadTime', 10e-6, ...
                              'B0', 2.98, ...
                              'adcSamplesLimit',8192);  
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