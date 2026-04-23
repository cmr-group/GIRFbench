function blip = get_blip(varargin)
    % Generate a triangle gradient blip waveform.
    %
    % Parameters (name-value pairs):
    %   dt       - raster time [s], default 10e-6
    %   n_slope  - number of ramp steps per side, default 10
    %   gmax     - max gradient amplitude [T/m], default 30e-3
    %   smax     - max slew rate [T/m/s], default 90
    %   sign     - polarity (+1 or -1), default 1

    parser = inputParser;
    addParameter(parser, 'dt', 10e-6);
    addParameter(parser, 'n_slope', 10);
    addParameter(parser, 'gmax', 30e-3);
    addParameter(parser, 'smax', 90);
    addParameter(parser, 'sign', 1);
    parse(parser, varargin{:});
    p = parser.Results;

    % Peak amplitude: limited by both gmax and slew rate
    a_peak = min(p.gmax, p.smax * p.n_slope * p.dt);
    if p.smax * p.n_slope * p.dt > p.gmax
        actual_slew = p.gmax / (p.n_slope * p.dt);
        warning('n_slope=%d exceeds gmax; slew derated to %.1f (%.1f%% of smax)', ...
            p.n_slope, actual_slew, 100 * actual_slew / p.smax);
    end

    ramp_up = linspace(0, a_peak, p.n_slope + 1);
    ramp_dn = linspace(a_peak, 0, p.n_slope + 1);
    blip = p.sign * [ramp_up, ramp_dn(2:end)];
end
