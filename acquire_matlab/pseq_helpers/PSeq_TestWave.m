classdef PSeq_TestWave < PSeq_Base
    % Base class for GIRF measurement test waves.

    % Supports thin slice and field camera measurements.  Derived classes should implement
    % prep_waves() to create the lists: all_test_waves, all_test_waves_neg, and
    % all_areas.

    properties
        do_adc;
        dt_adc;
        N_adc;
        total_duration;
        adc_delay;
        slew;

        adc_segments;
        adc_samples_per_segment;

        adc;
        duration_delay;

        N_waves;
        all_test_waves;
        all_test_waves_neg;
        all_areas;
    end

    methods
        function obj = PSeq_TestWave(pparams,varargin)
        % Construct.

        % Parameters
        % ----------
        % do_adc : bool, optional
        %     Should the ADC be played out, False for Skope measurements, True for thin
        %     slice, by default True
        % dt_adc : float, optional
        %     dwell time of ADC in [seconds], by default 4e-6
        % N_adc : int, optional
        %     number of points in ADC, by default 20000
        % total_duration : _type_, optional
        %     How long should this component be in [seconds], intended to get more
        %     consistent timings between differnet length test waves, by default 100e-3
        % adc_delay : float, optional
        %     Delay time for the ADC, in [seconds], by default 1e-3
        % slew : float, optional
        %     Slew rate override, in [T/m/s], if None keep using the current system setting,
        %     by default None

            obj = obj@PSeq_Base(pparams);

            p = inputParser;
            addParameter(p, 'do_adc', true);
            addParameter(p, 'dt_adc', 4e-6);
            addParameter(p, 'N_adc', 20000);
            addParameter(p, 'total_duration', 100e-3);
            addParameter(p, 'adc_delay', 1e-3);
            addParameter(p, 'slew', -1);

            parse(p, varargin{:});

            fields = fieldnames( p.Results );

            for n = 1:numel( fields )
	            obj.( fields{ n } ) = p.Results.( fields{ n } );
            end

            % -------------------------
            % Override slew rate if defined in argument
            if obj.slew > 0
                obj.slew = obj.pparams.sys.gamma * obj.slew;
            else
                obj.slew = obj.pparams.sys.maxSlew;
            end

            if obj.do_adc
                if obj.N_adc > obj.pparams.sys.adcSamplesLimit
                    [obj.adc_segments, obj.adc_samples_per_segment] = ...
                    mr.calcAdcSeg(obj.N_adc, obj.dt_adc, obj.pparams.sys);
                else
                    obj.adc_segments = 1;
                    obj.adc_samples_per_segment = obj.N_adc;
                end

                obj.adc = mr.makeAdc(obj.N_adc, 'dwell', obj.dt_adc, 'delay', obj.adc_delay, 'system', obj.pparams.sys);

                if mr.calcDuration(obj.adc) > obj.total_duration
                    disp('WARNING: ADC duration longer than total_duration, this is probably unintended.');                   
                end

            end

            obj.duration_delay = mr.makeDelay(obj.total_duration);
        end

        function prep_waves(obj)
            % Abstract base class for deriving the test waveforms.
            disp('WARNING: Base class prep_waves() should never be called.');  
        end


        function allblocks = build_blocks(obj, varargin)
            % Build all blocks to add to sequence.
            % 
            % Parameters
            % ----------
            % Parameters
            % ----------
            % idx : int, optional
            %     Index of the test waves to play.
            %     By default 0
            % polarity : int, optional
            %     Polarity of the test wave to play, either -1, 0, or 1.  0 does not play anything.
            %     By default 1
            % 
            % Returns
            % -------
            % cell array of cell arrays
            %     All blocks to add, outer list is for sequential blocks, inner blocks are all
            %     components to play within the block.

            parser = inputParser;
            addParameter(parser, 'idx', 1);
            addParameter(parser, 'polarity', 1.0);

            parse(parser, varargin{:});
            opt = parser.Results;

            grads_to_play = {};
            
            if opt.polarity == 1
                obj.all_test_waves{opt.idx}.channel = obj.pparams.channels{3};
                grads_to_play{end+1} = obj.all_test_waves{opt.idx};
            elseif opt.polarity == -1
                obj.all_test_waves_neg{opt.idx}.channel = obj.pparams.channels{3};
                grads_to_play{end+1} = obj.all_test_waves_neg{opt.idx};
            end

            if obj.do_adc
                if obj.pparams.rf_spoil
                    obj.adc.phaseOffset = obj.pparams.rf_spoil_phase;
                end
                grads_to_play{end+1} = obj.adc;
            end

            grads_to_play{end+1} = obj.duration_delay;

            allblocks = {grads_to_play};
            
        end

    end
end