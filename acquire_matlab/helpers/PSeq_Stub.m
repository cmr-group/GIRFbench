classdef PSeq_Stub < PSeq_Base
    % Sequence component that plays slice excitation and optional phase encoding
    % 
    % This is made for GIRF measurements, so the FOV is assumed to be square.

    properties

    end

    methods
        function obj = PSeq_Stub(pparams,varargin)

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

        end


        function allblocks = build_blocks(obj, idx1, idx2, offset)


        end

    end
end