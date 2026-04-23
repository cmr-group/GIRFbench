function bipolar = get_bipolar(varargin)
    % Generate a bipolar gradient waveform (triangle or trapezoid).
    % Ramps seamlessly from +peak through zero to -peak (no dwell at zero).
    % Zero net moment by symmetry.
    %
    % Shape (triangle):  0 → +a → -a → 0
    % Shape (trapezoid): 0 → +gmax ── +gmax → -gmax ── -gmax → 0
    %
    % Parameters (name-value pairs):
    %   dt       - raster time [s], default 10e-6
    %   t_lobe   - duration of one lobe (zero to zero) [s], default 1e-3
    %   gmax     - max gradient amplitude [T/m], default 30e-3
    %   smax     - max slew rate [T/m/s], default 90
    %   sign     - polarity of first lobe (+1 or -1), default 1

    parser = inputParser;
    addParameter(parser, 'dt', 10e-6);
    addParameter(parser, 't_lobe', 1e-3);
    addParameter(parser, 'gmax', 30e-3);
    addParameter(parser, 'smax', 90);
    addParameter(parser, 'sign', 1);
    parse(parser, varargin{:});
    p = parser.Results;

    % n_lobe = total samples for one lobe (zero to zero)
    n_lobe = max(2, round(p.t_lobe / p.dt));
    % n_half = half-lobe (used for triangle ramp)
    n_half = round(n_lobe / 2);
    % n_ramp = min samples to ramp from 0 to gmax at smax
    n_ramp = ceil(p.gmax / (p.smax * p.dt));

    if n_half <= n_ramp
        % Triangle at smax, peak below gmax
        a_peak = p.smax * n_half * p.dt;
        ramp_up = linspace(0, a_peak, n_half + 1);
        ramp_mid = linspace(a_peak, -a_peak, 2 * n_half + 1);
        ramp_dn = linspace(-a_peak, 0, n_half + 1);
        bipolar = p.sign * [ramp_up, ramp_mid(2:end), ramp_dn(2:end)];
    else
        % Long lobe: trapezoid at smax, peak = gmax
        n_flat = n_lobe - 2 * n_ramp;
        ramp_up = linspace(0, p.gmax, n_ramp + 1);
        flat_pos = p.gmax * ones(1, n_flat);
        ramp_mid = linspace(p.gmax, -p.gmax, 2 * n_ramp + 1);
        flat_neg = -p.gmax * ones(1, n_flat);
        ramp_dn = linspace(-p.gmax, 0, n_ramp + 1);
        bipolar = p.sign * [ramp_up, flat_pos, ramp_mid(2:end), flat_neg, ramp_dn(2:end)];
    end
end
