classdef Parameters
    % Class that describes all parameters to use in a GMLVQ run
    
    properties
        doztr (1,1) { GMLVQ.Helpers.mustBeLogical } = true
        mode (1,1) GMLVQ.Mode = GMLVQ.Mode.GMLVQNS
        rndinit (1,1) { GMLVQ.Helpers.mustBeLogical } = false
        mu (1,1) { mustBeNumeric } = 0
        decfac (1,1) { mustBeNumeric } = 1.5
        incfac (1,1) { mustBeNumeric } = 1.1
        ncop (1,1) { mustBeInteger, mustBeNonnegative } = 5
        etam (1,1) { mustBeNumeric } = NaN;
        etap (1,1) { mustBeNumeric } = NaN;
        rngseed (1,1) { mustBeNumeric } = 291024
        showlegend (1,1) { GMLVQ.Helpers.mustBeLogical } = true
        randomization (1,1) { mustBeNumeric } = 0.02
        useKMeans (1,1) { GMLVQ.Helpers.mustBeLogical } = true
        rocClass (1,1) { mustBeInteger, mustBeNonnegative } = 1
    end
    
    properties(Transient)
        warned = false;
    end
    
    methods
        % Set the parameters defined in key-value pairs in the
        % call-varargs and let them override the default values defined
        % above
        function params = Parameters(varargin)
            params = GMLVQ.Helpers.parseClassProperties(params, varargin{:});
            
            % Now do some validation and conditional parameter setting
            switch params.mode
                case {GMLVQ.Mode.GMLVQ, GMLVQ.Mode.GMLVQNS}
                    etam = 2;
                    etap = 1;
                    if params.mode == 0; disp('Matrix relevances without null-space correction'); end
                    if params.mode == 1; disp('Matrix relevances with null-space correction'); end
                case GMLVQ.Mode.GRLVQ
                    disp('Diagonal relevances, not encouraged, sensitive to step sizes');
                    etam = 0.2;
                    etap = 0.1;
                case GMLVQ.Mode.GLVQ
                    disp('GLVQ without relevances');
                    etam = 0;
                    etap = 1;
            end
            
            % Now check if etam and etap are set before applying defaults
            % RJV
            if isnan(params.etap) % Not set, apply default
                params.etap = etap;
            else
                if params.etap == etap
                    disp('Notice: setting etap to the default value has no effect');
                else
                    fprintf('Notice: etap set to %g instead of recommended %g for current mode\n', params.etap, etap);
                end
            end
            if isnan(params.etam) 
                params.etam = etam;
            else
                if params.etam == etam
                    disp('Notice: setting etam to the default value has no effect');
                else
                    fprintf('Notice: etam set to %g instead of recommended %g for current mode\n', params.etam, etam);
                end
            end
            
            params.sanityCheck();
        end
        
        function params = sanityCheck(params)
            if ~params.warned
                % Some general warnings
                if params.mode == GMLVQ.Mode.GLVQ && params.etam ~= 0
                    warning('When mode is set to GLVQ, etam should be set to 0, proceeding at your own risk.');
                end

                if ~params.doztr
                    disp('No z-score transformation, you may have to adjust step sizes');
                    if params.mode < 3; disp('Rescale relevances for proper interpretation'); end
                end            
                params.warned = true;
            else
                disp('...');
            end
        end
        
        function obj = unwarn(obj)
            obj.warned = false;
        end
        
       function obj = set.etam(obj,value)
            obj.etam = value;
            obj.unwarn();
%             disp('set-etam called');
        end
    end
end

