classdef PSeq_Excite_PE < PSeq_Base
    % Sequence component that plays slice excitation and optional phase encoding
    % 
    % This is made for GIRF measurements, so the FOV is assumed to be square.

    properties
        duration;
        thickness;
        tbw;
        flip;
        app;
        fov;
        N_pe;
        max_slew;
        do_refocus;
        do_prewind;

        slew;
        refocus_time;
        pe_areas;
        max_area;
        
        rfp;
        gss;
        gss_re;
        pe_grad1;
        pe_grad2;
        amp1;
        amp2;
    end

    methods
        function obj = PSeq_Excite_PE(pparams,varargin)

            obj = obj@PSeq_Base(pparams);

            p = inputParser;
            addParameter(p,'duration',6e-3);
            addParameter(p,'thickness',1e-3);
            addParameter(p,'tbw',4);
            addParameter(p,'flip',50);
            addParameter(p,'app',0.4);
            addParameter(p,'fov',200e-3);
            addParameter(p,'N_pe',0);
            addParameter(p,'max_slew',-1);
            addParameter(p,'do_refocus',true);
            addParameter(p,'do_prewind',true);

            parse(p, varargin{:});

            fields = fieldnames( p.Results );

            for n = 1:numel( fields )
	            obj.( fields{ n } ) = p.Results.( fields{ n } );
            end

            % ----------------------
            if obj.flip > pi
                obj.flip = obj.flip * pi / 180;
            end

            if obj.max_slew > 0
                obj.slew = obj.pparams.sys.gamma*maxSlew;
            else
                obj.slew = obj.pparams.sys.maxSlew;
            end


            % Get RF, slice select and refocus/prewinding gradient (assumes RF is coming out symmetric for prewinder)
            [obj.rfp, obj.gss, obj.gss_re] = mr.makeSincPulse(obj.flip, 'apodization', obj.app, 'duration', obj.duration, ...
                                                    'system', obj.pparams.sys, 'timeBwProduct', obj.tbw, 'delay', obj.pparams.sys.rfDeadTime, ...
                                                    'sliceThickness', obj.thickness, 'maxSlew', obj.slew, 'use', 'excitation');
            
            % Figure out the longest time needed between refocusing and phase encoding
            obj.refocus_time = 0;
            if obj.N_pe > 0
                obj.pe_areas = ((0:obj.N_pe-1) - floor(obj.N_pe/2))/obj.fov;
                obj.max_area = max(abs(obj.pe_areas));  % This gets used to rescale amplitude
                pe_temp = mr.makeTrapezoid(obj.pparams.channels{1}, 'area', obj.max_area, 'system', obj.pparams.sys);
                obj.refocus_time = mr.calcDuration(pe_temp);
            end
            if obj.do_refocus
                if mr.calcDuration(obj.gss_re) >= obj.refocus_time
                    obj.refocus_time = mr.calcDuration(obj.gss_re);
                else
                    % Remake gss_re to match longer phase encode time
                    obj.gss_re = mr.makeTrapezoid(obj.pparams.channels{3}, 'duration', obj.refocus_time, ...
                                                    'area', obj.gss_re.area, 'system', obj.pparams.sys);
                end
            end
               
            % Set up final phase encode gradients if needed     
            if obj.N_pe > 0
                obj.pe_grad1 = mr.makeTrapezoid(obj.pparams.channels{1}, 'duration', obj.refocus_time, ...
                                                        'area', obj.max_area, 'system', obj.pparams.sys);        
                
                obj.pe_grad2 = mr.makeTrapezoid(obj.pparams.channels{2}, 'duration', obj.refocus_time, ...
                                                        'area', obj.max_area, 'system', obj.pparams.sys);  
                
                obj.amp1 = obj.pe_grad1.amplitude;
                obj.amp2 = obj.pe_grad2.amplitude;
            end

        end

        function allblocks = build_blocks(obj, varargin)
            % Build all blocks to add to sequence.
            % 
            % Parameters
            % ----------
            % idx0, idx1 : int
            %     Index for the 0th channel area.  Will be ignored if the areas given for this
            %     channel were not inputted with multiple values.
            % 
            % offset: float
            %     Slice offset [m]
            % 
            % Returns
            % -------
            % cell array of cell arrays
            %     All blocks to add, outer list is for sequential blocks, inner blocks are all
            %     components to play within the block.

            parser = inputParser;
            addParameter(parser, 'idx1', 1);
            addParameter(parser, 'idx2', 1);
            addParameter(parser, 'seg_idx', -1);
            addParameter(parser, 'offset', 0);

            parse(parser, varargin{:});
            opt = parser.Results;

            obj.gss.channel = obj.pparams.channels{3};
            obj.rfp.freqOffset = obj.gss.amplitude * opt.offset;

            if obj.pparams.rf_spoil
                obj.rfp.phaseOffset = obj.pparams.rf_spoil_phase;
            end

            allblocks = {};

            % --- Slice select prewinder
            if obj.do_prewind
                obj.gss_re.channel = obj.pparams.channels{3};
                blocks = {obj.gss_re};
                if (opt.seg_idx >= 0)
                    blocks{end+1} = mr.makeLabel('SET', 'TRID', opt.seg_idx);
                end
                allblocks{end+1} = blocks;
            end

            % --- RF and slice select gradient
            blocks = {obj.rfp, obj.gss};
            if ((opt.seg_idx >= 0) && ~obj.do_prewind)
                blocks{end+1} = mr.makeLabel('SET', 'TRID', opt.seg_idx);
            end
            allblocks{end+1} = blocks;

            % --- Spatial encoding and slice select refocusing
            blocks = {};
            
            if obj.N_pe > 0
                area1 = obj.pe_areas(opt.idx1);
                area2 = obj.pe_areas(opt.idx2);
                
                obj.pe_grad1.channel = obj.pparams.channels{1};
                obj.pe_grad2.channel = obj.pparams.channels{2};
                
                obj.pe_grad1.amplitude = obj.amp1 * area1/obj.max_area;
                obj.pe_grad2.amplitude = obj.amp2 * area2/obj.max_area;
                
                blocks{end+1} = obj.pe_grad1;
                blocks{end+1} = obj.pe_grad2;
            end
                
                
            if obj.do_refocus
                obj.gss_re.channel = obj.pparams.channels{3};
                blocks{end+1} = obj.gss_re;
            end
                
            allblocks{end+1} = blocks;


        end

    end
end