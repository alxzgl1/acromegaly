classdef Helpers
    % Class with helper functions used in this library
    methods(Static)
        
        % Function to override the default values specified with key-value pairs
        function options = parseParameters(defaults, varargin)
            optionNames = fieldnames(defaults); %Fix RJV change fieldnames to properties
            nArgs = length(varargin);
            if mod(nArgs, 2) ~= 0
                error('Set key/value pairs as input parameters')
            end
            
            for pair = reshape(varargin, 2, [])
                name = pair{1};
                
                if any(strcmp(name, optionNames))
                    defaults.(name) = pair{2};
                else
                    error('%s is not a valid parameter name', name)
                end
            end
            
            options = defaults;
        end
        
        % Function to wrap the default-override from parseParameters for classes
        function obj = parseClassProperties(obj, varargin)
            % First create a struct of the fields, to be used in the other function
            warning('off', 'MATLAB:structOnObject');
            defaults = struct(obj);
            names = fieldnames(defaults);
            warning('on', 'MATLAB:structOnObject');
            
            % Remove empty fields (without default values, they will be calculated later)
            defaults = rmfield(defaults, names(structfun(@isempty, defaults))); %RJV prevented setting etam and etap
            
            % Call other function to parse the parameters
            params = GMLVQ.Helpers.parseParameters(defaults, varargin{:});
            
            % Put the parameters in the class object
            names = fieldnames(params);
            for i = 1:numel(names)
                obj.(names{i}) = params.(names{i});
            end
        end
        
        % Function to set RNG correctly
        function setRNG(rngseed)
            rng('default');
            rng(rngseed);
        end
        
        % Function to check for being logical
        function mustBeLogical(a)
            if ~islogical(a)
                error('Value must be logical.');
            end
        end
        
        % Function to calculate the ROC
        function [tpr, fpr, auroc, thresholds] = roc(binlbl, score, nThresholds)

            % Modification April 2016: 

            % threshold-based computation of the ROC
            % scores are rescaled to fall into the range 0...1 
            % and then compared with nthresh equi-distant thresholds
            % note that nthresh must be large enough to guarantee correct
            % ROC computation 
            % planned improvement: 
            %  - re-scaling according to observed range of scores
            %  - efficient rank-based computation of ROC, mapped to thresholds 
            %    in order to facilitate threshold-averages

            if nargin < 3 || isempty(nThresholds)
                nThresholds = 5000;
            end

            % Heuristic rescaling of the scores to range 0...1
            % TODO: Re-scale according to observed range values

            score = 1 ./ (1 + exp(score/2));

            target = binlbl';
            length(binlbl);

            tu = unique(target);
            t1 = tu(1);
            t2 = tu(2);

            % For proper "threshold-averages" ROC (see paper by Fawcett) we use "nthresh" equi-distant
            % thresholds between 0 and 1

            if length(binlbl) > 1250; nThresholds = 4*length(binlbl); end
            if mod(nThresholds, 2) == 1; nThresholds = nThresholds - 1; end % Make sure it is even
            thresholds = linspace(0, 1, nThresholds + 1);

            fpr = zeros(1, nThresholds + 1);
            tpr = fpr;
            tpr(1) = 1; fpr(1) = 1;

            for i = 1:nThresholds - 1
                % Count true positives, false positives
                tp = sum(target(score > thresholds(i+1)) == t2);
                fp = sum(target(score > thresholds(i+1)) == t1);
                fn = sum(target(score <= thresholds(i+1)) == t2);
                tn = sum(target(score <= thresholds(i+1)) == t1);

                % Compute corresponding rates
                tpr(i+1) = tp / (tp + fn);
                fpr(i+1) = fp / (tn + fp);
            end

            % Simple numerical integration
            auroc = -trapz(fpr, tpr); % Minus sign due to order of fpr, tpr

        end
    end
end

