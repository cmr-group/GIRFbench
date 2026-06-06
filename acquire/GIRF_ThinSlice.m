%%
clear all
clc

BASE_DIR = 'D:\Dropbox\projects\matlab2\';

addpath([BASE_DIR 'pulseq151\matlab\']);
addpath('helpers\');

%%
% Here we set some general parameters for the sequence

TR = 500e-3;
N_av = 4;

smax0 = 110;  % T/m/s
gmax0 = 40;  % mT/m

FOVz = 100e-3;
Nslices = 5;

params = PSeq_Params('max_grad', gmax0, 'max_slew', smax0);
pseq0 = PSeq_Base(params);

gmax = 0.99*params.sys.maxGrad;
smax = 0.99*params.sys.maxSlew;
dt = params.sys.gradRasterTime;

% Slice offsets
slice_shift = 2e-3;  % Helps to have some "signal" in all slices for debugging
all_offsets = linspace(-FOVz/2, FOVz/2, Nslices) + slice_shift;

%%
% All waves are just arbitrary gradient waveforms in Hz/m

all_waves = {};

% ------- Blips
all_waves{end+1} = get_blip('dt', dt, 'n_slope', 4,...
                'gmax', gmax, 'smax', smax);

all_waves{end+1} = get_blip('dt', dt, 'n_slope', 18,...
                'gmax', gmax, 'smax', smax);

all_waves{end+1} = get_blip('dt', dt, 'n_slope', 26,...
                'gmax', gmax, 'smax', smax);

% ------- Chirps
% Lowest frequency, low freq at end, 10ms
all_waves{end+1} = get_chirp('f1', 0, 'f2', 5000, 't_chirp', 10e-3, 'dt', dt, ...
                    'reverse', true, 'gmax', gmax, 'smax', smax, 'sweep_exp', 2.0, ...
                    'max_k', 500, 'auto_margin', true, 'zero_moment', false);

% Lowest frequency, low freq at start, 20ms
all_waves{end+1} = get_chirp('f1', 0, 'f2', 20000, 't_chirp', 20e-3, 'dt', dt, ...
                    'reverse', false, 'gmax', gmax, 'smax', smax, 'sweep_exp', 1.0, ...
                    'max_k', 500, 'auto_margin', true, 'zero_moment', false);

% Lowest frequency, low freq at end, refocused at end, 10ms
all_waves{end+1} = get_chirp('f1', 0, 'f2', 10000, 't_chirp', 10e-3, 'dt', dt, ...
                    'reverse', true, 'gmax', gmax, 'smax', smax, 'sweep_exp', 1.0, ...
                    'max_k', 500, 'auto_margin', true, 'zero_moment', true);

% ------- Bipolars
all_waves{end+1} = get_bipolar('dt', dt, 't_lobe', 0.5e-3,...
                    'gmax', gmax, 'smax', smax);

all_waves{end+1} = get_bipolar('dt', dt, 't_lobe', 1.0e-3,...
                    'gmax', gmax, 'smax', smax);

%%
% PSeq_GIRFGeneral is a helper class that handles almost all of the GIRF measurement, check the function
% docstring to see all of the options that are supported

pseq_girf = PSeq_GIRFGeneral(params, 'test_waves_in', all_waves, 'excite_mode', 'Thin');

%%

pseq0.reinit_seq();
iTR = 0;

for i_av = 1:N_av
for channels = {{'y', 'z', 'x'}, {'x', 'z', 'y'}, {'x', 'y', 'z'}}

for wave_idx = 1:numel(all_waves)
for polarity = [-1, 1]
    
    params.channels = channels{1};
    
    pseq0.track_time = 0;
    pseq0.add_block_list( pseq_girf.build_blocks('wave_idx', wave_idx, 'wave_polarity', polarity) );

    req_delay = TR-pseq0.track_time;
    if req_delay > 0
        pseq0.add_delay(req_delay);
    else
        disp('WARNING: The current TR was longer than the prescribed TR')
    end
    iTR = iTR + 1;

end  % polarity
end  % wave_idx

end  % channels
end  % av


%%
% check whether the timing of the sequence is correct
% ----------
[ok, error_report]=pseq0.seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% ----------
seq_duration = seconds(pseq0.seq.duration);
seq_duration.Format = 'hh:mm:ss';
fprintf('Sequence Duration: %s \n', seq_duration)


%%
% Write out the sequence and a bunch of metadata with it
% ----------

% --- Set name
dtime = datetime(); dtime.Format = 'yyyy_MMdd';
seq_name = sprintf('GIRF_Thin_%s',  dtime);

% Increment the filename so we never overwrite a sequence
i_name = 0;
out_fname = sprintf('../acquire/export/%s_%03d.seq', seq_name, i_name);
while isfile(out_fname)
    i_name = i_name + 1;
    out_fname = sprintf('../acquire/export/%s_%03d.seq', seq_name, i_name);
end
seq_name = sprintf('%s_%03d', seq_name, i_name);
disp(out_fname)

% --- Save params to .json
out_params.N_av = N_av;
out_params.TR = TR;
out_params.N_TR = iTR;
out_params.N_Waves = numel(all_waves);
out_params.smax = smax0;
out_params.gmax = gmax0;

out_params.all_offsets = all_offsets;
out_params.FOVz = FOVz;
out_params.N_slices = Nslices;

out_params.seq_name = sprintf('%s.seq', seq_name);
out_params.seq_duration = sprintf('%s', seq_duration);
out_params = orderfields(out_params);

out_params.general_params = struct();
out_params.general_params = pseq_girf.get_scan_info(out_params.general_params);

jsonStr = jsonencode(out_params, 'PrettyPrint', true);
writelines(jsonStr, sprintf('../acquire/export/%s.json', seq_name));


% --- Save sequence
pseq0.seq.setDefinition('FOV', [60e-3, 60e-3, 60e-3]);
pseq0.seq.setDefinition('Name', seq_name);

pseq0.seq.write(out_fname);   % Output sequence for scanner
save(sprintf('../acquire/export/%s_waves.mat', seq_name), 'all_waves');

fprintf('Done!\n');
