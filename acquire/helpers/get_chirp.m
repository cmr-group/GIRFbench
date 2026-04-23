function chirp = get_chirp(varargin)
    % Generate a chirp gradient waveform that sweeps from f1 to f2,
    % constrained by maximum amplitude (gmax) and slew rate (smax).
    %
    % Parameters (name-value pairs):
    %   dt      - raster time [s], default 10e-6
    %   t_chirp - duration [s], default 30e-3
    %   f1      - start frequency [Hz], default 0
    %   f2      - end frequency [Hz], default 15000
    %   gmax          - max gradient amplitude [T/m], default 30e-3
    %   smax          - max slew rate [T/m/s], default 90
    %   slew_margin   - safety factor for slew envelope (0-1), default 0.95
    %   smooth_dur    - envelope smoothing window duration [s], default 1e-3
    %   reverse       - flip the waveform (high freq first), default false
    %   sweep_exp     - frequency sweep exponent (1=linear, >1=more time at low freq), default 1
    %   max_k         - max k-space excursion [1/m], disabled when <= 0, default -1
    %   k_margin      - safety factor for k-space envelope (0-1), default 0.9
    %   auto_margin   - iteratively optimize margins to hit constraints, default false
    %   zero_moment   - null the zeroth gradient moment (append rewinder), default false
    %   rewinder_pos  - where to place the rewinder: 'end' or 'start', default 'end'

    parser = inputParser;
    addParameter(parser, 'dt', 10e-6);
    addParameter(parser, 't_chirp', 30e-3);
    addParameter(parser, 'f1', 0);
    addParameter(parser, 'f2', 15000);
    addParameter(parser, 'gmax', 30e-3);
    addParameter(parser, 'smax', 90);
    addParameter(parser, 'slew_margin', 0.95);
    addParameter(parser, 'smooth_dur', 1e-3);
    addParameter(parser, 'reverse', false);
    addParameter(parser, 'sweep_exp', 1);
    addParameter(parser, 'max_k', -1);
    addParameter(parser, 'k_margin', 0.9);
    addParameter(parser, 'auto_margin', false);
    addParameter(parser, 'zero_moment', false);
    addParameter(parser, 'rewinder_pos', 'end');
    parse(parser, varargin{:});
    p = parser.Results;

    GAMMA = 42.577478461e6; % Hz/T

    tt = 0:p.dt:p.t_chirp;
    k = p.sweep_exp;
    T = p.t_chirp;

    % Power-law frequency sweep: f(t) = f1 + (f2-f1)*(t/T)^k
    % k=1 is linear, k>1 spends more time at low frequencies
    f_inst = p.f1 + (p.f2 - p.f1) * (tt / T).^k;
    phase = 2 * pi * (p.f1 * tt + (p.f2 - p.f1) * T / (k + 1) * (tt / T).^(k + 1));

    % df/dt for the product-rule slew correction
    df_dt = k * (p.f2 - p.f1) / T * (tt / T).^(k - 1);
    alpha = 1 ./ sqrt(1 + (df_dt ./ (2 * pi * f_inst.^2 + eps)).^2);

    % Smoothing kernel (precompute)
    n_smooth = max(3, round(p.smooth_dur / p.dt));
    kern = hann(n_smooth);
    kern = kern / sum(kern);

    % Auto-margin: start at 1.0 and iteratively correct
    if p.auto_margin
        slew_margin = 1.0;
        k_margin = 1.0;
        n_iter = 5;
    else
        slew_margin = p.slew_margin;
        k_margin = p.k_margin;
        n_iter = 1;
    end

    for iter = 1:n_iter
        % Amplitude envelope with slew constraint
        amp_env = min(p.gmax, slew_margin * alpha .* p.smax ./ (2 * pi * f_inst + eps));

        % K-space constraint
        if p.max_k > 0
            amp_env = min(amp_env, k_margin * p.max_k * 2 * pi * f_inst);
        end

        % Smooth the envelope transition
        amp_env = conv(amp_env, kern, 'same');
        amp_env = min(amp_env, p.gmax);

        chirp = amp_env .* sin(phase);

        % Truncate at the last zero crossing so the waveform ends near zero
        sign_changes = find(chirp(1:end-1) .* chirp(2:end) <= 0);
        if ~isempty(sign_changes)
            last = sign_changes(end);
            if abs(chirp(last)) <= abs(chirp(last+1))
                chirp = chirp(1:last);
            else
                chirp = chirp(1:last+1);
            end
        end
        chirp(end+1) = 0;

        % Measure actual constraints and update margins
        if p.auto_margin && iter < n_iter
            if p.max_k > 0
                if p.reverse
                    actual_k = max(abs(cumsum(fliplr(chirp)) * p.dt));
                else
                    actual_k = max(abs(cumsum(chirp) * p.dt));
                end
                k_margin = k_margin * (p.max_k / actual_k);
                k_margin = min(k_margin, 1); % ensure k margin doesn't exceed 1
            end

            actual_slew = max(abs(diff(chirp) / p.dt));
            % fprintf('iter: %d  actual_slew: %.2f\n', iter, actual_slew/GAMMA)
            slew_margin = 0.999999 * slew_margin * (p.smax / actual_slew);
        end
    end

    if p.reverse
        chirp = fliplr(chirp);
    end

    % Append rewinder to null residual moment (always at the end)
    if p.zero_moment
        M0 = sum(chirp) * p.dt;
        if abs(M0) > eps
            area = abs(M0);
            s = -sign(M0);
            % Try triangle (minimum duration): area = smax * t_ramp^2
            t_ramp = sqrt(area / p.smax);
            a_peak = p.smax * t_ramp;
            if a_peak > p.gmax
                % Need trapezoid: ramp to gmax, hold, ramp down
                t_ramp = p.gmax / p.smax;
                t_flat = area / p.gmax - t_ramp;
                n_ramp = max(1, ceil(t_ramp / p.dt));
                n_flat = max(0, round(t_flat / p.dt));
                ramp_up = linspace(0, s*p.gmax, n_ramp+1);
                flat = s * p.gmax * ones(1, n_flat);
                ramp_dn = linspace(s*p.gmax, 0, n_ramp+1);
                rewinder = [ramp_up, flat, ramp_dn(2:end)];
            else
                n_ramp = max(1, ceil(t_ramp / p.dt));
                ramp_up = linspace(0, s*a_peak, n_ramp+1);
                ramp_dn = linspace(s*a_peak, 0, n_ramp+1);
                rewinder = [ramp_up, ramp_dn(2:end)];
            end
            % Scale rewinder to get exact zero moment
            rew_moment = sum(rewinder) * p.dt;
            if abs(rew_moment) > eps
                rewinder = rewinder * (-M0 / rew_moment);
            end
            if strcmp(p.rewinder_pos, 'start')
                chirp = [rewinder, chirp];
            else
                chirp = [chirp, rewinder];
            end
        end
    end

    % Verify constraints
    slew = diff(chirp) / p.dt;
    fprintf('Max slew: %.2f  Smax: %.2f\n', max(abs(slew))/GAMMA, p.smax/GAMMA);
    if max(abs(slew)) > p.smax
        warning('Chirp slew rate (%.1f) exceeds smax (%.1f)', max(abs(slew))/GAMMA, p.smax/GAMMA);
        figure()
        plot(abs(slew)/GAMMA);
    end
    if p.max_k > 0
        actual_k = max(abs(cumsum(chirp) * p.dt));
        fprintf('Max k: %.4f  max_k: %.4f\n', actual_k, p.max_k);
    end
    if p.zero_moment
        fprintf('Final moment (M0): %.2e\n', sum(chirp) * p.dt);
    end
    if p.auto_margin
        fprintf('Auto margins: slew_margin=%.4f  k_margin=%.4f\n', slew_margin, k_margin);
    end
end
