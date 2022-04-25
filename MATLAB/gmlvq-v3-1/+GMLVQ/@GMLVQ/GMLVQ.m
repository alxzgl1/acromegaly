classdef GMLVQ < matlab.mixin.Copyable
    % Base GMLVQ class which is the starting point of all things to do with this library.
    
    properties
        data GMLVQ.DataPair
        totalsteps (1,1) { mustBeInteger }
        plbl (:,1) { mustBeInteger } %RJV TODO: Make dependent public set/get method with private var
        params (1,1) GMLVQ.Parameters
    end
    
    properties (SetAccess = private) %FIX RJV Do not save/load these as they are stored with plbl
        nPrototypes = 0
        nClasses = 0
    end
    
    properties (Dependent)
        nDimensions
        nFeatureVectors
    end
    
    methods
        
        % Constructor
        % @param featureVectors NxM matrix of scalars
        % @param labels Nx1 vector of integers
        % @param parameters GMLVQ.Parameters optional
        % @param totalsteps integer optional, default = 10
        % @param prototypeLabels vector of integers optional
        % @out GMLVQ.GMLVQ
        function obj = GMLVQ(featureVectors, labels, parameters, totalsteps, prototypeLabels)
            obj.data = GMLVQ.DataPair(featureVectors, labels);
            % Check data consistency
            GMLVQ.DataPair.mustBeIncreasing(obj.data.labels);
            
            % Handle empty parameters
            if nargin < 3 || isempty(parameters) || ~isa(parameters, 'GMLVQ.Parameters')
                obj.params = GMLVQ.Parameters();
                warning('Defaulted to default parameters!');
            else
                obj.params = parameters;
            end
            
            % Handle empty totalsteps
            if nargin < 4 || isempty(totalsteps)
                obj.totalsteps = 10;
                disp('Defaulted to 10 training steps');
            else
                obj.totalsteps = totalsteps;
            end
            
            % Handle empty plbl
            if nargin < 5 || isempty(prototypeLabels)
                obj.plbl = 1:length(unique(obj.data.labels));
                disp('Defaulted to one prototype per class');
            else
                obj.plbl = prototypeLabels;
            end
            disp('Prototype configuration:'); disp(obj.plbl');
            
            % Check whether we have a consistent model
            obj.checkConsistency();
        end
        
        % Setter function for the plbl that also updates its dependencies
        function set.plbl(gmlvq, v)
            gmlvq.plbl = v;
            gmlvq.nPrototypes = length(v);
            gmlvq.nClasses = length(unique(v));
        end
        
        % Getter for nDimensions
        function res = get.nDimensions(this)
            res = this.data.nDimensions;
        end
        
        % Getter for nFeatureVectors
        function res = get.nFeatureVectors(this)
            res = this.data.nFeatureVectors;
        end
        
        % These functions are in other files
        gmlvq_result = runSingle(gmlvq);
        res = runValidation(this, nRuns, percentage);
        res = runL1O(this)
    end
    
    methods (Access = protected)
        % These functions are in other files
        checkConsistency(gmlvq);
    end
end

