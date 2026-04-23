classdef PSeq_Base  < handle
    % Base class and common functions for pyPulseq helper classes.

    properties
       pparams;
       seq;
       track_time = 0.0;
       track_total_time = 0.0;
    end

    methods
        function obj = PSeq_Base(pparams)
            % PSeq_Base Construct an instance of this class
      
            obj.pparams = pparams;
            obj.seq = mr.Sequence(obj.pparams.sys);
        end

        function reinit_seq(obj)
            % Reset the sequence.

            obj.pparams.rf_spoil_idx = 0;
            obj.pparams.rf_spoil_phase = 0;
            obj.track_time = 0;
            obj.track_total_time = 0;
            obj.seq = mr.Sequence(obj.pparams.sys);
        end

        function add_delay(obj, delay)
            % Add a simple delay to the sequence.
            % 
            % Parameters
            % ----------
            % delay : float
            %     Time of the delay in [seconds]

            obj.seq.addBlock(mr.makeDelay(delay));
            obj.track_time = obj.track_time + delay; 
            obj.track_total_time = obj.track_total_time + delay; 
        end

        function out = get_seq_time(obj)
            % Return the time in seconds of the entire sequence.
            
            out = obj.seq.duration();
        end

        function add_block_list(obj, all_blocks)
            % Add a list of block elements to this sequence.
            % 
            % See each components build_blocks function for argument list.

            for i = 1:numel(all_blocks)
                % Add block to sequence
                obj.seq.addBlock(all_blocks{i});
            end

            dur = duration_block_list(all_blocks);
            obj.track_time = obj.track_time + dur;
            obj.track_total_time = obj.track_total_time + dur; 
        end

       
    end
end

function dur = duration_block_list(all_blocks)
    dur = 0;
    for i = 1:numel(all_blocks)
        dur = dur + mr.calcDuration(all_blocks{i});
    end
end