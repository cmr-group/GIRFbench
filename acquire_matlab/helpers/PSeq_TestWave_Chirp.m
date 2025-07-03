classdef PSeq_TestWave_Chirp < PSeq_TestWave
    % Base class for GIRF measurement test waves.

    % Supports thin slice and field camera measurements.  Derived classes should implement
    % prep_waves() to create the lists: all_test_waves, all_test_waves_neg, and
    % all_areas.

    properties
        all_chirps;
        wave_delay;
    end

    methods
        function obj = PSeq_TestWave_Chirp(pparams,varargin)
        % Construct.

        % Parameters
        % ----------
        % wave_delay : float, optional
        %     Delay time [seconds] before playing out the waveform, to allow some sampling
        %     of the ADC before any gradients are played out, by default 2e-3

        % Notes
        % -----
        % A lot of this is hard coded for initial testing, this needs to be much more
        % customizeable and a final design strategy selected after testing how the different
        % chirps perform.

            obj = obj@PSeq_TestWave(pparams, varargin{:});

            p = inputParser;
            p.KeepUnmatched=true;  % Allows for passthrough inhereted options to base class

            addParameter(p, 'wave_delay', 2e-3);

            parse(p, varargin{:});

            fields = fieldnames( p.Results );

            for n = 1:numel( fields )
	            obj.( fields{ n } ) = p.Results.( fields{ n } );
            end

            obj.prep_waves()
        end


        function prep_waves(obj)
            
            obj.all_chirps = cell(1,1);
            obj.all_chirps{1} = get_chirp('dt', obj.pparams.sys.gradRasterTime, ...
                                          'gmax', obj.pparams.sys.maxGrad, ...
                                          'smax', obj.slew, ...
                                          'max_k', -1, ...
                                          'median_k', -1);

            obj.N_waves = numel(obj.all_chirps);

            for i = 1:obj.N_waves
                wave = mr.makeArbitraryGrad( ...
                    obj.pparams.channels{3}, ...
                    obj.all_chirps{i}, ...
                    'delay', obj.wave_delay, ...
                    'first', 0, ...
                    'last', 0, ...
                    'system', obj.pparams.sys);

                obj.all_test_waves{i} = wave;
                obj.all_areas(i) = wave.area;

                wave = mr.makeArbitraryGrad( ...
                    obj.pparams.channels{3}, ...
                    -obj.all_chirps{i}, ...
                    'delay', obj.wave_delay, ...
                    'first', 0, ...
                    'last', 0, ...
                    'system', obj.pparams.sys);

                obj.all_test_waves_neg{i} = wave;

            end

        end


    end  % methods
end % classdef


function chirp = get_chirp(varargin)
    % Design a chirp waveform subject to moment constraints.

    % Parameters
    % ----------
    % dt : float
    %     Raster time [seconds]
    % t_chirp : float
    %     Duration of the chirp [seconds]
    % f1 : float
    %     Starting frequency [Hz]
    % f2 : float
    %     Ending frequency [Hz]
    % max_k : float, optional
    %     Maximum absolute M0 [1/m]
    % median_k : int, optional
    %     Maximum median M0 [1/m]
    % gmax : float, optional
    %     Maximum gradient strength [T/m], by default 30e-3
    % smax : int, optional
    %     Maximum slew rate [T/m/s], by default 90

    % Returns
    % -------
    % array
    %     The chirp waveform [T/m] with raster time dt
    
    parser = inputParser;
    addParameter(parser, 't_chirp', 30e-3);
    addParameter(parser, 'dt', 10e-6);
    addParameter(parser, 'f1', 0);
    addParameter(parser, 'f2', 15000);
    addParameter(parser, 'gmax', 30e-3);
    addParameter(parser, 'smax', 90);
    addParameter(parser, 'max_k', -1);
    addParameter(parser, 'median_k', -1);


    parse(parser, varargin{:});
    opt = parser.Results;

    tt = 0:opt.dt:opt.t_chirp;
    ft = opt.f1 + (opt.f2 - opt.f1) .* tt ./ opt.t_chirp;

    gmax_vec = ones(1, numel(tt)) * opt.gmax;

    Gct = gmax_vec .* sin(2 .* pi .* (opt.f1 .* tt + (opt.f2 - opt.f1) .* tt.^2 ./ 2 ./ opt.t_chirp));
    senv = 0.98 .* opt.smax ./ (2 .* pi .* gmax_vec .* ft + 1e-64);
    senv(senv>1) = 1;
    
    chirp = senv .* Gct;

    [chirp_m0_min, chirp_m0_max, chirp_m0_median] = get_m0(chirp, opt.dt);
    chirp_ok = (max(abs([chirp_m0_min, chirp_m0_max])) < opt.max_k || opt.max_k<0) && (chirp_m0_median < opt.median_k || opt.median_k<0);
    % fprintf('OK: %d  min_m0 = %.1f  max_m0 = %.1f  median_m0 = %.1f \n', chirp_ok, chirp_m0_min, chirp_m0_max, chirp_m0_median);

    idx_first_lobe = find(sin(2 .* pi .* (opt.f1 .* tt + (opt.f2 - opt.f1) .* tt.^2 ./ 2 ./ opt.t_chirp)) < 0, 1);

    lobe_reduce = 1.0;
    while ~chirp_ok
        lobe_reduce = lobe_reduce * 0.95;
        
        gmax_vec = ones(1, numel(tt)) * opt.gmax;
        gmax_vec(1:idx_first_lobe) = lobe_reduce*opt.gmax;

        k = kaiser(round(0.9*idx_first_lobe), 9);
        k = k./sum(k);

        gmax_vec = [ones(1,numel(k))*gmax_vec(1), gmax_vec, ones(1,numel(k))*gmax_vec(end)];
        gmax_vec = conv(gmax_vec, k, 'full');
        gmax_vec = gmax_vec(ceil(1.5*numel(k)):end);
        gmax_vec = gmax_vec(1:end-numel(gmax_vec)+numel(tt));

        
        Gct = gmax_vec .* sin(2 .* pi .* (opt.f1 .* tt + (opt.f2 - opt.f1) .* tt.^2 ./ 2 ./ opt.t_chirp));
        senv = 0.98 .* opt.smax ./ (2 .* pi .* gmax_vec .* ft + 1e-64);
        senv(senv>1) = 1;
        
        chirp = senv .* Gct;
        
        [chirp_m0_min, chirp_m0_max, chirp_m0_median] = get_m0(chirp, opt.dt);
        chirp_ok = (max(abs([chirp_m0_min, chirp_m0_max])) < opt.max_k || opt.max_k<0) && (chirp_m0_median < opt.median_k || opt.median_k<0);
        % fprintf('OK: %d  min_m0 = %.1f  max_m0 = %.1f  median_m0 = %.1f \n', chirp_ok, chirp_m0_min, chirp_m0_max, chirp_m0_median);

    end

    slew_check = diff(chirp)./opt.dt;
    if max(abs(slew_check)) > opt.smax
        disp('ERROR: Chirp generation gave slew rate > smax.')
    end
end

function [min_m0, max_m0, median_m0] = get_m0(chirp, dt)
    m0 = dt * cumsum(chirp);
    min_m0 = min(m0);
    max_m0 = max(m0);
    median_m0 = median(m0);
end