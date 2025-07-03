classdef PSeq_Moments3D < PSeq_Base
    % Sequence component that plays traps on each axis (or fewer).
    % 
    % This is made to be fairly flexible in designing the shortest possible time gradients,
    % to capture a given starting moment or ending moment.  The intended use is at the end
    % of a TR for refocusing/spoiling, but it can also be used for phase encoding or other uses.

    properties
        start_areas;
        end_areas;
        start_polarities;
        end_polarities;

        max_areas;
        active_grad;

        max_duration;

        all_grads;
        all_amps;
        all_modes;
    end

    methods
        function obj = PSeq_Moments3D(pparams,varargin)
        % Construct.
        % 
        % Parameters
        % ----------
        % start_areas : length 3 list-like
        %     Areas [1/m] before the start of this component for each gradient channel.  Can be
        %     None, a single float, or a list-like of floats.  If None, it assumes there
        %     will be no gradient played on that axis.
        %     By default (None, None, None)
        % end_areas : length 3 list-like
        %     Areas [1/m] desired at the end of this component for each gradient channel.  Uses
        %     the same format as start_areas.
        %     By default (None, None, None)
        % start_polarities : length 3 list-like
        %     Lists possible polarities for each start areas.  Mostly for GIRF measurements,
        %     otherwise this is just 1
        % end_polarities : length 3 list-like
        %     Lists possible polarities for each end areas.  Mostly for GIRF measurements, 
        %     otherwise this is just 1
        % 
        % Notes
        % -----
        % If the areas provided for starting or ending are a list, the calls to run this
        % component will be indexed, using the given list of areas.  Currently this class is
        % planned to only handle a single list for each direction, but this can be expanded
        % later if needed (TODO).
        % 
        % All trapezoids are created with the same duration, matching the minimum duration
        % needed for largest required area.  For phase-encoding lists of areas, the
        % amplitude only is scaled, timing remains the same (TODO: Add options for shortest
        % durations on each channel?)

            obj = obj@PSeq_Base(pparams);

            p = inputParser;
            addParameter(p,'start_areas', {[], [], []});
            addParameter(p,'end_areas', {[], [], []});
            addParameter(p,'start_polarities', {1, 1, 1});
            addParameter(p,'end_polarities', {1, 1, 1});

            parse(p, varargin{:});

            fields = fieldnames( p.Results );

            for n = 1:numel( fields )
	            obj.( fields{ n } ) = p.Results.( fields{ n } );
            end

            % -------------------------
            obj.max_areas = zeros(1,3);
            for i = 1:3
                obj.max_areas(i) = get_max_difference(obj.start_areas{i}, obj.end_areas{i}, obj.start_polarities{i}, obj.end_polarities{i});
            end

            if ~(any(obj.max_areas))
                disp('WARNING: PSeq_Moments3D found no areas needed in any direction, will play .1ms delay instead')
            end
    
            obj.active_grad = obj.max_areas > 0;

            % Get minimum time gradient for each channel, to get final duration
            temp_grads = {};
            for i = 1:3
                if obj.active_grad(i)
                    grad_ = mr.makeTrapezoid( ...
                        obj.pparams.channels{i}, ...
                        'area', obj.max_areas(i), ...
                        'system', obj.pparams.sys);
                    temp_grads{end+1} = grad_;
                end
            end
            
            if any(obj.max_areas)
                obj.max_duration = mr.calcDuration(temp_grads);
                
                obj.all_grads = cell(1,3);
                obj.all_amps = zeros(1,3);
                for i = 1:3
                    if obj.active_grad(i)
                        grad_ = mr.makeTrapezoid( ...
                                obj.pparams.channels{i}, ...
                                'duration', obj.max_duration, ...
                                'area', obj.max_areas(i), ...
                                'system', obj.pparams.sys);

                        obj.all_grads{i} = grad_;
                        obj.all_amps(i) = grad_.amplitude;
                    else
                        obj.all_grads{i} = [];
                        obj.all_amps(i) = 0;
                    end
                end


                obj.all_modes = cell(1,3);
                for i = 1:3
                    if obj.active_grad(i)
                        if numel(obj.start_areas{i}) > 1 && numel(obj.end_areas{i}) > 1
                            disp('ERROR: PSeq_Moments3D does not currently handle lists for both start and end areas');
                            obj.all_modes{i} = 'none';
                        end

                        if isscalar(obj.start_areas{i}) && isscalar(obj.end_areas{i})
                            obj.all_modes{i} = 'scalar';
                        elseif numel(obj.start_areas{i}) > 1
                            obj.all_modes{i} = 'start';
                        elseif numel(obj.end_areas{i}) > 1
                            obj.all_modes{i} = 'end';
                        end
                    else
                        obj.all_modes{i} = 'none';
                    end
                end

            end
        end


        function allblocks = build_blocks(obj, varargin)
            % Build all blocks to add to sequence.
            % 
            % Parameters
            % ----------
            % idx0, idx1, idx2 : int
            %     Index for the 0th channel area.  Will be ignored if the areas given for this
            %     channel were not inputted with multiple values.
            % polarities : floats
            %     Polarities for start and end areas for each channel.  This has been added for
            %     GIRF measurements, allows the areas to be reversed.
            % 
            % Returns
            % -------
            % cell array of cell arrays
            %     All blocks to add, outer list is for sequential blocks, inner blocks are all
            %     components to play within the block.

            parser = inputParser;
            addParameter(parser, 'idx1', 1);
            addParameter(parser, 'idx2', 1);
            addParameter(parser, 'idx3', 1);
            addParameter(parser, 'start_pol1', 1.0);
            addParameter(parser, 'start_pol2', 1.0);
            addParameter(parser, 'start_pol3', 1.0);
            addParameter(parser, 'end_pol1', 1.0);
            addParameter(parser, 'end_pol2', 1.0);
            addParameter(parser, 'end_pol3', 1.0);

            parse(parser, varargin{:});
            opt = parser.Results;

            all_idx = [opt.idx1, opt.idx2, opt.idx3];
            all_start_pol = [opt.start_pol1, opt.start_pol2, opt.start_pol3];
            all_end_pol = [opt.end_pol1, opt.end_pol2, opt.end_pol3];
            
            grads_to_play = {};

            for i = 1:3
                if obj.active_grad(i)
                    if strcmp(obj.all_modes{i}, 'start')
                        start_area = all_start_pol(i) * obj.start_areas{i}(all_idx(i));
                    else
                        start_area = all_start_pol(i) * obj.start_areas{i};
                    end

                    if strcmp(obj.all_modes{i}, 'end')
                        end_area = all_end_pol(i) * obj.end_areas{i}(all_idx(i));
                    else
                        end_area = all_end_pol(i) * obj.end_areas{i};
                    end

                    obj.all_grads{i}.amplitude = obj.all_amps(i) * (end_area - start_area) / obj.max_areas(i);
                    obj.all_grads{i}.channel = obj.pparams.channels{i};

                    grads_to_play{end+1} = obj.all_grads{i};
                end
            end

            % This is a (time-wasting) hack when no gradient is needed, TODO: handle better.
            if numel(grads_to_play) == 0
                grads_to_play = {mr.makeDelat(.1e-3)};
            end
                
            allblocks = {grads_to_play};
        end

    end
end

function max_diff = get_max_difference(start_moment, end_moment, start_pol, end_pol)
    % Get the maximum difference (end-start).
    % 
    % start_moment and end_moment may be [], a scalar, or a list/array, and the maximum difference of
    % all combinations is returned.  If either argument is [], will return 0
    
    if (numel(start_moment) == 0) || (numel(end_moment) == 0) || (numel(start_pol) == 0) || (numel(end_pol) == 0)
        max_diff = 0;
        return
    end

    max_diff = 0;

    for i0 = 1:numel(start_moment)
    for i1 = 1:numel(end_moment)
    for i2 = 1:numel(start_pol)
    for i3 = 1:numel(end_pol)
        diff = end_pol(i3) * end_moment(i1) - start_pol(i2) * start_moment(i0);
        max_diff = max([abs(diff), max_diff]);
    end
    end
    end
    end

end