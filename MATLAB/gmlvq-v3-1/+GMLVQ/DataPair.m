classdef DataPair < matlab.mixin.Copyable
    % Class that represents a set of feature vectors with their corresponding labels
    
    properties
        featureVectors (:,:) { mustBeNumeric, mustBeFinite, mustBeReal }
        labels (:,1) { mustBeInteger } = []
    end
    
    properties (SetAccess = private)
        nDimensions = 0
        nFeatureVectors = 0
    end
    
    methods
        % Constructor
        % @param featureVectors NxM matrix
        % @param labels Nx1 vector with corresponding labels
        % @out GMLVQ.DataPair
        function dataSet = DataPair(featureVectors, labels)
            dataSet.featureVectors = featureVectors;
            dataSet.labels = labels;
            
            if(size(dataSet.labels, 2) > 1)
                dataSet.labels = dataSet.labels';
            end
        end
        
        % Setter function for the data that also updates its dependencies
        function set.featureVectors(dataSet, fv)
            dataSet.featureVectors = fv;
            dataSet.nDimensions = size(fv, 2);
            dataSet.nFeatureVectors = size(fv, 1);
        end
    end
    
    methods (Static)
        
        % Performs the label value check (only use 1, 2, ..., n as class names)
        function mustBeIncreasing(x)
            if ~isempty(x) && (min(x) ~= 1 || max(x) ~= length(unique(x)))
                error('Data labels should be 1, 2, 3, ..., n!');
            end
        end
        
    end
end