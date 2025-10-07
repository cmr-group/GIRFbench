classdef PSeq_TestWave_Blip < PSeq_TestWave
    % Blip test waveforms for GIRF measurement

    % Supports thin slice and field camera measurements.  Derived classes should implement
    % prep_waves() to create the lists: all_test_waves, all_test_waves_neg, and
    % all_areas.

    properties
        all_chirps;
        wave_delay;
        min_ramp;
        ramp_inc;
    end

    methods
        function obj = PSeq_TestWave_Blip(pparams,varargin)
        % Construct.

        % Parameters
        % ----------
        % min_ramp : float, optional
        %     Ramp time [seconds] for the smallest blip, by default 60e-6
        % ramp_inc : float, optional
        %     Time [seconds] to increase each subsequent ramp time in N_waves, by default 10e-6
        % wave_delay : float, optional
        %     Delay time [seconds] before playing out the waveform, to allow some sampling
        %     of the ADC before any gradients are played out, by default 2e-3
        % N_waves : int, optional
        %     Nuber of waveforms to generate, by default 12

            obj = obj@PSeq_TestWave(pparams, varargin{:});

            p = inputParser;
            p.KeepUnmatched=true;  % Allows for passthrough inhereted options to base class

            addParameter(p, 'min_ramp', 60e-6);
            addParameter(p, 'ramp_inc', 10e-6);
            addParameter(p, 'N_waves', 12);
            addParameter(p, 'wave_delay', 2e-3);

            parse(p, varargin{:});

            fields = fieldnames( p.Results );

            for n = 1:numel( fields )
	            obj.( fields{ n } ) = p.Results.( fields{ n } );
            end
            
            obj.min_ramp = ceil(obj.min_ramp/obj.pparams.sys.gradRasterTime)*obj.pparams.sys.gradRasterTime;
            obj.ramp_inc = ceil(obj.ramp_inc/obj.pparams.sys.gradRasterTime)*obj.pparams.sys.gradRasterTime;

            obj.prep_waves()
        end


        function prep_waves(obj)

            for i = 1:obj.N_waves

                ramp_time = obj.min_ramp + (i-1) * obj.ramp_inc;

                % .99999 due to float rounding errors sometimes exceeding slew rate limit
                amp = 0.99999 * obj.slew * ramp_time;

                wave = mr.makeTrapezoid( ...
                    obj.pparams.channels{3}, ...
                    'flatTime', 0, ...
                    'amplitude', amp, ...
                    'riseTime', ramp_time, ...
                    'fallTime', ramp_time, ...
                    'delay', obj.wave_delay, ...
                    'system', obj.pparams.sys);

                obj.all_test_waves{i} = wave;
                obj.all_areas(i) = wave.area;

                wave = mr.makeTrapezoid( ...
                    obj.pparams.channels{3}, ...
                    'flatTime', 0, ...
                    'amplitude', -amp, ...
                    'riseTime', ramp_time, ...
                    'fallTime', ramp_time, ...
                    'delay', obj.wave_delay, ...
                    'system', obj.pparams.sys);

                obj.all_test_waves_neg{i} = wave;

            end

        end


    end  % methods
end % classdef