classdef PSeq_GIRFGeneral < PSeq_Base

    properties
        % Excitation options
        excite_mode
        rf_duration
        rf_thickness
        rf_tbw
        rf_flip
        rf_app
        do_refocus
        do_prewind

        % PE options
        PE_fov
        N_pe

        % ADC options
        do_adc
        dt_adc
        N_adc

        % Timing options
        acq_order
        excite_delay
        adc_delay
        wave_delay
        spoil_delay

        % Fixed timing options
        fixedtime_excite
        fixedtime_adc
        fixedtime_wave
        fixedtime_spoiler

        % Spoiling options
        spoil_mode
        spoil_moment

        % General options
        min_meas_time
        test_waves_in
        max_slew

        % Derived / computed properties
        slew
        rfp
        rfp_delay
        gss
        gss_re
        trig
        trig_delay
        adc
        adc_segments
        adc_samples_per_segment
        N_waves
        all_test_waves
        all_areas
        all_durations
        all_test_waves_neg
        max_test_time
        pe_areas
        pe_max_area
        pe_grad1
        pe_grad2
        pe_amp1
        pe_amp2
        spoil_max_area
        spoil_grad1
        spoil_grad2
        spoil_grad3
        spoil_amp1
        spoil_amp2
        spoil_amp3
    end

    methods
        function obj = PSeq_GIRFGeneral(pparams, options)
            % PSeq_GIRFGeneral  Construct a GIRF-general pulse sequence object.
            %
            %   obj = PSeq_GIRFGeneral(pparams)
            %   obj = PSeq_GIRFGeneral(pparams, 'flip', 90, 'N_pe', 64)
            %
            %   Inputs:
            %       pparams     - Parameter struct passed to PSeq_Base
            %
            %   Name-Value Arguments:
            %
            %     Excitation:
            %       excite_mode  - Excitation mode: 'Skope', 'Thin', 'Thin_PE', or 'NonSel' (default: 'Thin')
            %       rf_duration  - RF pulse duration in seconds (default: 6e-3)
            %       rf_thickness - Slice thickness in meters (default: 1e-3)
            %       rf_tbw       - Time-bandwidth product (default: 4)
            %       rf_flip      - Flip angle in degrees or radians (default: 50)
            %       rf_app       - Apodization factor (default: 0.4)
            %       do_refocus   - Apply slice-select refocusing (default: true)
            %       do_prewind   - Apply slice-select prewinding (default: true)
            %
            %     Phase encoding:   (only used when excite_mode == 'Thin_PE')
            %       PE_fov       - Phase-encode field of view in meters (default: 240e-3)
            %       N_pe         - Number of phase encodes, 0 to disable (default: 0)
            %
            %     ADC:
            %       do_adc       - Enable ADC readout (default: true)
            %       dt_adc       - ADC dwell time in seconds (default: 4e-6)
            %       N_adc        - Number of ADC samples (default: 20000)
            %
            %     Timing:
            %       acq_order    - 'excite_before' or 'excite_after' (default: 'excite_before')
            %       excite_delay - Delay before excitation in seconds (default: 0)
            %       adc_delay    - Delay before ADC in seconds (default: 0)
            %       wave_delay   - Delay before test waveform in seconds (default: 0.2e-3)
            %       spoil_delay  - Delay before spoiler in seconds (default: 0)
            %                        (Note that this is *after* "min_meas_time" which is its own form of spoiler delay)
            %
            %     Fixed timing (overrides other timing, -1 to disable):
            %       fixedtime_excite  - Fixed excitation time (default: -1)
            %       fixedtime_adc     - Fixed ADC time (default: -1)
            %       fixedtime_wave    - Fixed waveform time (default: -1)
            %       fixedtime_spoiler - Fixed spoiler time (default: -1)
            %
            %     Spoiling:
            %       spoil_mode   - 'none', 'refocus', or 'spoil' (default: 'refocus')
            %       spoil_moment - Spoiler gradient moment (default: 3000)
            %
            %     General:
            %       min_meas_time - Minimum measurement time in seconds (default: 110e-3)
            %       test_waves_in - Cell array of test waveforms (default: [])
            %       max_slew      - Max slew rate in T/m/s, -1 for system default (default: -1)
            
            arguments
                pparams

                % Excitation options
                options.excite_mode {mustBeMember(options.excite_mode, {'Skope', 'Thin', 'Thin_PE', 'NonSel'})} = 'Thin'
                options.rf_duration = 6e-3
                options.rf_thickness = 1e-3
                options.rf_tbw = 4
                options.rf_flip = 50
                options.rf_app = 0.4
                options.do_refocus = true
                options.do_prewind = true

                % PE options
                options.PE_fov = 240e-3
                options.N_pe = 0

                % ADC options
                options.do_adc = true
                options.dt_adc = 4e-6
                options.N_adc = 20000

                % Timing options
                options.acq_order {mustBeMember(options.acq_order, {'excite_before', 'excite_after'})} = 'excite_before'
                options.excite_delay = 0
                options.adc_delay = 0
                options.wave_delay = 0.2e-3
                options.spoil_delay = 0e-3
                
                % Options for fixed timing (overrides other timing settings)
                options.fixedtime_excite = -1
                options.fixedtime_adc = -1
                options.fixedtime_wave = -1
                options.fixedtime_spoiler = -1

                % Spoiling options
                options.spoil_mode {mustBeMember(options.spoil_mode, {'none', 'refocus', 'spoil'})} = 'refocus'
                options.spoil_moment = 3000

                % General options
                options.min_meas_time = 110e-3
                options.test_waves_in = []
                options.max_slew = -1

            end

            obj = obj@PSeq_Base(pparams);

            fields = fieldnames(options);
            for n = 1:numel(fields)
                obj.(fields{n}) = options.(fields{n});
            end

            if obj.max_slew > 0
                obj.slew = obj.pparams.sys.gamma*obj.max_slew;
            else
                obj.slew = obj.pparams.sys.maxSlew;
            end

            % -------- Prep Excitation --------
            if obj.excite_mode == "Skope"
                obj.prep_excite_skope();
            elseif obj.excite_mode == "Thin"
                obj.prep_excite_thin();
            elseif obj.excite_mode == "Thin_PE"
                obj.prep_excite_thin_pe();
            elseif obj.excite_mode == "NonSel"
                obj.prep_excite_nonselective();
            end

            % -------- Prep ADC --------
            if obj.do_adc
                if obj.N_adc > obj.pparams.sys.adcSamplesLimit
                    [obj.adc_segments, obj.adc_samples_per_segment] = ...
                    mr.calcAdcSeg(obj.N_adc, obj.dt_adc, obj.pparams.sys);
                else
                    obj.adc_segments = 1;
                    obj.adc_samples_per_segment = obj.N_adc;
                end

                obj.adc = mr.makeAdc(obj.N_adc, 'dwell', obj.dt_adc, 'delay', obj.adc_delay, 'system', obj.pparams.sys);
            end

            % -------- Prep Test Waves --------
            obj.N_waves = numel(obj.test_waves_in);
            for i = 1:obj.N_waves
                wave = mr.makeArbitraryGrad( ...
                    obj.pparams.channels{3}, ...
                    obj.test_waves_in{i}, ...
                    'delay', obj.wave_delay, ...
                    'first', 0, ...
                    'last', 0, ...
                    'system', obj.pparams.sys);

                obj.all_test_waves{i} = wave;
                obj.all_areas(i) = wave.area;
                obj.all_durations(i) = mr.calcDuration(wave);

                wave = mr.makeArbitraryGrad( ...
                    obj.pparams.channels{3}, ...
                    -obj.test_waves_in{i}, ...
                    'delay', obj.wave_delay, ...
                    'first', 0, ...
                    'last', 0, ...
                    'system', obj.pparams.sys);

                obj.all_test_waves_neg{i} = wave;
            end
            obj.max_test_time = max(obj.all_durations);

            % -------- Prep Spoilers --------
            obj.spoil_max_area = max(abs(obj.all_areas));
            if ~isempty(obj.pe_areas)
                obj.spoil_max_area = max(obj.spoil_max_area, obj.pe_max_area);
            end
            
            if obj.spoil_mode == "spoil"
                obj.spoil_max_area = obj.spoil_max_area + obj.spoil_moment;
            end

            fprintf('Max spoiler area needed: %f \n', obj.spoil_max_area);

            obj.spoil_grad1 = mr.makeTrapezoid(obj.pparams.channels{1}, 'area', obj.spoil_max_area, 'maxSlew', obj.slew, 'system', obj.pparams.sys);
            obj.spoil_grad2 = mr.makeTrapezoid(obj.pparams.channels{2}, 'area', obj.spoil_max_area, 'maxSlew', obj.slew, 'system', obj.pparams.sys);
            obj.spoil_grad3 = mr.makeTrapezoid(obj.pparams.channels{3}, 'area', obj.spoil_max_area, 'maxSlew', obj.slew, 'system', obj.pparams.sys);

            obj.spoil_amp1 = obj.spoil_grad1.amplitude;
            obj.spoil_amp2 = obj.spoil_grad2.amplitude;
            obj.spoil_amp3 = obj.spoil_grad3.amplitude;
        
        end

        function prep_excite_skope(obj)
            obj.trig = mr.makeDigitalOutputPulse('ext1', 'duration', obj.pparams.sys.gradRasterTime);
            obj.trig_delay = mr.makeDelay(200e-6);
        end

        function prep_excite_nonselective(obj)
            if obj.rf_flip > pi
                obj.rf_flip = obj.rf_flip * pi / 180;
            end

            obj.rfp = mr.makeBlockPulse(obj.rf_flip, 'Duration', 0.3e-3, 'system', obj.pparams.sys, ...
                                        'delay', obj.pparams.sys.rfDeadTime, 'use','excitation');
            obj.rfp_delay = mr.makeDelay(0.5e-3);
        end

        function prep_excite_thin(obj)
            if obj.rf_flip > pi
                obj.rf_flip = obj.rf_flip * pi / 180;
            end

            % Get RF, slice select and refocus/prewinding gradient (assumes RF is coming out symmetric for prewinder)
            [obj.rfp, obj.gss, obj.gss_re] = mr.makeSincPulse(obj.rf_flip, 'apodization', obj.rf_app, 'duration', obj.rf_duration, ...
                                                    'system', obj.pparams.sys, 'timeBwProduct', obj.rf_tbw, 'delay', obj.pparams.sys.rfDeadTime, ...
                                                    'sliceThickness', obj.rf_thickness, 'maxSlew', obj.slew, 'use', 'excitation');
        end

        function prep_excite_thin_pe(obj)
            if obj.rf_flip > pi
                obj.rf_flip = obj.rf_flip * pi / 180;
            end

            % Get RF, slice select and refocus/prewinding gradient (assumes RF is coming out symmetric for prewinder)
            [obj.rfp, obj.gss, obj.gss_re] = mr.makeSincPulse(obj.rf_flip, 'apodization', obj.rf_app, 'duration', obj.rf_duration, ...
                                                    'system', obj.pparams.sys, 'timeBwProduct', obj.rf_tbw, 'delay', obj.pparams.sys.rfDeadTime, ...
                                                    'sliceThickness', obj.rf_thickness, 'maxSlew', obj.slew, 'use', 'excitation');


            obj.pe_areas = ((0:obj.N_pe-1) - floor(obj.N_pe/2))/obj.PE_fov;
            obj.pe_max_area = max(abs(obj.pe_areas));
            pe_temp = mr.makeTrapezoid(obj.pparams.channels{1}, 'area', obj.pe_max_area, 'maxSlew', obj.slew, 'system', obj.pparams.sys);
            refocus_time = mr.calcDuration(pe_temp);

            if obj.do_refocus
                if mr.calcDuration(obj.gss_re) >= refocus_time
                    refocus_time = mr.calcDuration(obj.gss_re);
                else
                    % Remake gss_re to match longer phase encode time
                    obj.gss_re = mr.makeTrapezoid(obj.pparams.channels{3}, 'duration', refocus_time, ...
                                                    'area', obj.gss_re.area, 'maxSlew', obj.slew, 'system', obj.pparams.sys);
                end
            end

            obj.pe_grad1 = mr.makeTrapezoid(obj.pparams.channels{1}, 'duration', refocus_time, ...
                                                            'area', obj.pe_max_area, 'maxSlew', obj.slew, 'system', obj.pparams.sys);        
                    
            obj.pe_grad2 = mr.makeTrapezoid(obj.pparams.channels{2}, 'duration', refocus_time, ...
                                                    'area', obj.pe_max_area, 'maxSlew', obj.slew, 'system', obj.pparams.sys);  
            
            obj.pe_amp1 = obj.pe_grad1.amplitude;
            obj.pe_amp2 = obj.pe_grad2.amplitude;
        end

        function allblocks = add_excite(obj, allblocks, options)
            if obj.excite_mode == "Skope"
                allblocks{end+1} = {obj.trig, obj.trig_delay};
            elseif obj.excite_mode == "NonSel"
                if obj.pparams.rf_spoil
                    obj.rfp.phaseOffset = obj.pparams.rf_spoil_phase;
                end
                allblocks{end+1} = {obj.rfp, obj.rfp_delay};
            elseif obj.excite_mode == "Thin" || obj.excite_mode == "Thin_PE"
                if obj.pparams.rf_spoil
                    obj.rfp.phaseOffset = obj.pparams.rf_spoil_phase;
                end
                obj.gss.channel = obj.pparams.channels{3};
                obj.rfp.freqOffset = obj.gss.amplitude * options.offset;

                % --- RF Prewinder if enabled
                if obj.do_prewind
                    obj.gss_re.channel = obj.pparams.channels{3};
                    allblocks{end+1} = {obj.gss_re};
                end

                % --- RF and slice select
                allblocks{end+1} = {obj.rfp, obj.gss};

                % --- Spatial encoding and slice select refocusing
                blocks = {};
                
                if obj.N_pe > 0 && obj.excite_mode == "Thin_PE"
                    area1 = obj.pe_areas(options.pe_idx1);
                    area2 = obj.pe_areas(options.pe_idx2);
                    
                    obj.pe_grad1.channel = obj.pparams.channels{1};
                    obj.pe_grad2.channel = obj.pparams.channels{2};
                    
                    obj.pe_grad1.amplitude = obj.pe_amp1 * area1/obj.pe_max_area;
                    obj.pe_grad2.amplitude = obj.pe_amp2 * area2/obj.pe_max_area;
                    
                    blocks{end+1} = obj.pe_grad1;
                    blocks{end+1} = obj.pe_grad2;   
                end
                    
                if obj.do_refocus
                    obj.gss_re.channel = obj.pparams.channels{3};
                    blocks{end+1} = obj.gss_re;
                end
                
                if numel(blocks) > 0
                    allblocks{end+1} = blocks;
                end

            end


        end

        function allblocks = build_blocks(obj, options)
            % build_blocks  Build all blocks to add to sequence.
            %
            %   allblocks = obj.build_blocks()
            %   allblocks = obj.build_blocks('pe_idx1', 3, 'offset', 0.01)
            %
            %   Name-Value Arguments:
            %       pe_idx1    - Index for the 1st PE channel area (default: 1)
            %       pe_idx2    - Index for the 2nd PE channel area (default: 1)
            %       offset  - Slice offset in meters (default: 0)
            %       wave_idx - Index for the waveform (default: 1)
            %       wave_polarity - Polarity of the waveform (default: 1)
            %
            %   Returns:
            %       allblocks - Cell array of cell arrays; outer list is sequential
            %                   blocks, inner cells are components within each block.
            %
            %  TODO: This currently doesnt handle the fixed timing options, which may require restructuring how blocks are built.  
            %         For now, if fixed timing is desired, build_blocks should be overridden in a child class to enforce the timing.
            %         That will make it so long this should be broken up
            arguments
                obj
                options.pe_idx1 = 1
                options.pe_idx2 = 1
                options.offset = 0

                options.wave_idx = 1
                options.wave_polarity = 1
            end

            allblocks = {};

            if obj.acq_order == "excite_before"
                % Add excitation block(s) 
                % =========================================
                allblocks = obj.add_excite(allblocks, options);

                % This builds up ADC and test wave
                % =========================================
                grads_to_play = {};

                % ADC
                % -------------------------
                if obj.do_adc
                    if obj.pparams.rf_spoil
                        obj.adc.phaseOffset = obj.pparams.rf_spoil_phase;
                    end
                    grads_to_play{end+1} = obj.adc;
                end
                
                % Get test wave if being played, but dont add it yet, as it may need to be combined with the spoiler
                % -------------------------
                if options.wave_idx > 0  % wave_idx = 0 means dont play a test waveform
                    if options.wave_polarity == 1
                        wave_to_play = obj.all_test_waves{options.wave_idx};
                    elseif options.wave_polarity == -1
                        wave_to_play = obj.all_test_waves_neg{options.wave_idx};
                    else
                        wave_to_play = [];
                    end
                else
                    wave_to_play = [];
                end
                
                if ~isempty(wave_to_play)
                    wave_to_play.channel = obj.pparams.channels{3};
                    wave_to_play.delay = obj.wave_delay;  % This should already be set in init
                    area_to_refocus = wave_to_play.area;
                    grads_to_play{end+1} = wave_to_play;
                else
                    area_to_refocus = 0;
                end

                grads_to_play{end+1} = mr.makeDelay(obj.min_meas_time);  % This ensures consistent timing even if no wave is played
                allblocks{end+1} = grads_to_play;

                
                % Spoilers/Refocusers
                % =========================================
                grads_to_play = {};

                if obj.spoil_mode == "spoil" || obj.spoil_mode == "refocus"
                    
                    % --- Set up all spoiler/refocusing gradients
                    if obj.excite_mode == "Thin_PE"
                        refocus_area1 = -obj.pe_areas(options.pe_idx1);
                        refocus_area2 = -obj.pe_areas(options.pe_idx2);
                    else
                        refocus_area1 = 0;
                        refocus_area2 = 0;
                    end
                    refocus_area3 = -area_to_refocus;

                    if obj.spoil_mode == "spoil"
                        refocus_area1 = refocus_area1 + obj.spoil_moment;
                        refocus_area2 = refocus_area2 + obj.spoil_moment;   
                        refocus_area3 = refocus_area3 + obj.spoil_moment;
                    end

                    obj.spoil_grad1.amplitude = obj.spoil_amp1 * refocus_area1/obj.spoil_max_area;
                    obj.spoil_grad2.amplitude = obj.spoil_amp2 * refocus_area2/obj.spoil_max_area;
                    obj.spoil_grad3.amplitude = obj.spoil_amp3 * refocus_area3/obj.spoil_max_area;

                    obj.spoil_grad1.channel = obj.pparams.channels{1};
                    obj.spoil_grad2.channel = obj.pparams.channels{2};
                    obj.spoil_grad3.channel = obj.pparams.channels{3};

                    obj.spoil_grad1.delay = obj.spoil_delay;
                    obj.spoil_grad2.delay = obj.spoil_delay;
                    obj.spoil_grad3.delay = obj.spoil_delay;

                    % --- Spoiler in x and y are straightforward: play if nonzero
                    if obj.spoil_grad1.amplitude ~= 0
                        grads_to_play{end+1} = obj.spoil_grad1;
                    end
                    
                    if obj.spoil_grad2.amplitude ~= 0
                        grads_to_play{end+1} = obj.spoil_grad2;
                    end

                    if obj.spoil_grad3.amplitude ~= 0
                        grads_to_play{end+1} = obj.spoil_grad3;
                    end

                end

                if ~isempty(grads_to_play)
                    allblocks{end+1} = grads_to_play;
                end
                

            elseif obj.acq_order == "excite_after"
                % Add test wave alone
                % Note: we keep the duration here constant to the longest test wave to ensure consistent timing, even if no wave is being played.
                % =========================================
                if options.wave_idx > 0  % wave_idx = 0 means dont play a test waveform
                    if options.wave_polarity == 1
                        wave_to_play = obj.all_test_waves{options.wave_idx};
                    elseif options.wave_polarity == -1
                        wave_to_play = obj.all_test_waves_neg{options.wave_idx};
                    else
                        wave_to_play = [];
                    end
                else
                    wave_to_play = [];
                end
                
                if ~isempty(wave_to_play)
                    wave_to_play.channel = obj.pparams.channels{3};
                    wave_to_play.delay = obj.wave_delay;  % This should already be set in init
                    area_to_refocus = wave_to_play.area;
                    % allblocks{end+1} = {wave_to_play, mr.makeDelay(obj.max_test_time)};
                    allblocks{end+1} = {wave_to_play};
                else
                    area_to_refocus = 0;
                    % allblocks{end+1} = {mr.makeDelay(obj.max_test_time)};
                    allblocks{end+1} = {mr.makeDelay(1e-3)};
                end


                % Add excitation block(s) 
                % =========================================
                allblocks = obj.add_excite(allblocks, options);

            
                % Add ADC 
                % =========================================
                grads_to_play = {};

                if obj.do_adc
                    if obj.pparams.rf_spoil
                        obj.adc.phaseOffset = obj.pparams.rf_spoil_phase;
                    end
                    grads_to_play{end+1} = obj.adc;
                end

                grads_to_play{end+1} = mr.makeDelay(obj.min_meas_time);  % This ensures consistent timing even if no wave is played
                allblocks{end+1} = grads_to_play;

                % Add Spoilers/Refocusers
                % =========================================
                grads_to_play = {};

                if obj.spoil_mode == "spoil" || obj.spoil_mode == "refocus"
                    if obj.excite_mode == "Thin_PE"
                        refocus_area1 = -obj.pe_areas(options.pe_idx1);
                        refocus_area2 = -obj.pe_areas(options.pe_idx2);
                    else
                        refocus_area1 = 0;
                        refocus_area2 = 0;
                    end
                    refocus_area3 = -area_to_refocus;

                    if obj.spoil_mode == "spoil"
                        refocus_area1 = refocus_area1 + obj.spoil_moment;
                        refocus_area2 = refocus_area2 + obj.spoil_moment;   
                        refocus_area3 = refocus_area3 + obj.spoil_moment;
                    end

                    obj.spoil_grad1.amplitude = obj.spoil_amp1 * refocus_area1/obj.spoil_max_area;
                    obj.spoil_grad2.amplitude = obj.spoil_amp2 * refocus_area2/obj.spoil_max_area;
                    obj.spoil_grad3.amplitude = obj.spoil_amp3 * refocus_area3/obj.spoil_max_area;

                    obj.spoil_grad1.channel = obj.pparams.channels{1};
                    obj.spoil_grad2.channel = obj.pparams.channels{2};
                    obj.spoil_grad3.channel = obj.pparams.channels{3};

                    obj.spoil_grad1.delay = obj.spoil_delay;
                    obj.spoil_grad2.delay = obj.spoil_delay;
                    obj.spoil_grad3.delay = obj.spoil_delay;

                    if obj.spoil_grad1.amplitude ~= 0
                        grads_to_play{end+1} = obj.spoil_grad1;
                    end
                    
                    if obj.spoil_grad2.amplitude ~= 0
                        grads_to_play{end+1} = obj.spoil_grad2;
                    end

                    if obj.spoil_grad3.amplitude ~= 0
                        grads_to_play{end+1} = obj.spoil_grad3;
                    end
                end

                if ~isempty(grads_to_play)
                    allblocks{end+1} = grads_to_play;
                end
   

            end

        end

        function s = get_scan_info(obj, s)
            % get_scan_info  Add scan description fields to a struct.
            %
            %   s = obj.get_scan_info(s)
            %
            %   Populates s with scalar/string fields describing the scan
            %   configuration. Intended for JSON export — no large arrays.

            % Excitation options
            s.excite_mode = obj.excite_mode;
            s.rf_duration = obj.rf_duration;
            s.rf_thickness = obj.rf_thickness;
            s.rf_tbw = obj.rf_tbw;
            s.rf_flip_rad = obj.rf_flip;
            s.rf_flip_deg = obj.rf_flip * 180 / pi;
            s.rf_apodization = obj.rf_app;
            s.do_refocus = obj.do_refocus;
            s.do_prewind = obj.do_prewind;

            % Phase encoding
            s.PE_fov_m = obj.PE_fov;
            s.N_pe = obj.N_pe;
            if obj.N_pe > 0 && ~isempty(obj.pe_max_area)
                s.pe_max_area = obj.pe_max_area;
            end

            % ADC
            s.do_adc = obj.do_adc;
            s.dt_adc_s = obj.dt_adc;
            s.N_adc = obj.N_adc;
            if obj.do_adc
                s.adc_segments = obj.adc_segments;
                s.adc_samples_per_segment = obj.adc_samples_per_segment;
                s.adc_duration_s = obj.N_adc * obj.dt_adc;
            end

            % Timing
            s.acq_order = obj.acq_order;
            s.excite_delay_s = obj.excite_delay;
            s.adc_delay_s = obj.adc_delay;
            s.wave_delay_s = obj.wave_delay;
            s.spoil_delay_s = obj.spoil_delay;
            s.min_meas_time_s = obj.min_meas_time;

            % Fixed timing
            s.fixedtime_excite = obj.fixedtime_excite;
            s.fixedtime_adc = obj.fixedtime_adc;
            s.fixedtime_wave = obj.fixedtime_wave;
            s.fixedtime_spoiler = obj.fixedtime_spoiler;

            % Spoiling
            s.spoil_mode = obj.spoil_mode;
            s.spoil_moment = obj.spoil_moment;
            s.spoil_max_area = obj.spoil_max_area;

            % Test waveforms (summary, not the arrays)
            s.N_waves = obj.N_waves;
            s.max_test_time_s = obj.max_test_time;
            if ~isempty(obj.all_areas)
                s.wave_area_min = min(obj.all_areas);
                s.wave_area_max = max(obj.all_areas);
            end
            if ~isempty(obj.all_durations)
                s.wave_duration_min_s = min(obj.all_durations);
                s.wave_duration_max_s = max(obj.all_durations);
            end

            % Slew rate
            s.max_slew_setting = obj.max_slew;
            s.effective_slew = obj.slew;

            % System parameters from pparams
            s.sys_maxGrad = obj.pparams.sys.maxGrad;
            s.sys_maxSlew = obj.pparams.sys.maxSlew;
            s.sys_gradRasterTime = obj.pparams.sys.gradRasterTime;
            s.sys_rfRasterTime = obj.pparams.sys.rfRasterTime;
            s.sys_rfDeadTime = obj.pparams.sys.rfDeadTime;
            s.sys_adcSamplesLimit = obj.pparams.sys.adcSamplesLimit;
            s.rf_spoil = obj.pparams.rf_spoil;
        end

    end  % methods
end  % classdef
