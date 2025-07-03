%%
addpath('D:/Dropbox/projects/matlab2/pulseq150/matlab/');
addpath('./helpers/');

%%
clearvars
clc


%%
specs.MaxGrad = 32;
specs.MaxSlew = 110;

TR = 500e-3;
N_av = 4;

params = PSeq_Params(specs);

pseq0 = PSeq_Base(params);

pseq_excite = PSeq_Excite_PE(params);

pseq_test = PSeq_TestWave_Chirp(params, 'slew', 100);

pseq_refocus = PSeq_Moments3D(params, ...
    'start_areas', {[], [], pseq_test.all_areas}, ...
    'end_areas', {0, 0, 0}, ...
    'start_polarities', {1, 1, [-1, 0, 1]} );

%%
% Slice offsets
FOVz = 80e-3;
Nslices = 5;
slice_shift = 2e-3;  % Helps to have some "signal" in all slices for debugging
all_offsets = linspace(-FOVz/2, FOVz/2, Nslices) + slice_shift;

%% Make Seq
% =================================================
%
pseq0.reinit_seq();

for i_av = 1:N_av
for channels = {{'y', 'z', 'x'}, {'x', 'z', 'y'}, {'x', 'y', 'z'}}
for polarity = [-1, 0, 1]
for offset = all_offsets
 
    params.channels = channels{1};
    params.increment_rf_spoiling();
    
    pseq0.track_time = 0;
    pseq0.add_block_list( pseq_excite.build_blocks('offset', offset) );
    pseq0.add_block_list( pseq_test.build_blocks('polarity', polarity) );
    pseq0.add_block_list( pseq_refocus.build_blocks('end_pol3', polarity) );

    req_delay = TR-pseq0.track_time;
    if req_delay > 0
        pseq0.add_delay(req_delay);
    end

end 
end
end
end

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

% ----------
pseq0.seq.setDefinition('FOV', [320e-3 320e-3 FOVz]);
pseq0.seq.setDefinition('Name', 'temp_GIRF');
pseq0.seq.setDefinition('MaxAdcSegmentLength', pseq_test.adc_samples_per_segment);

pseq0.seq.write('export/phantom_girf.seq');   % Output sequence for scanner
fprintf('Done!\n');