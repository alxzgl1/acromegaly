classdef (Abstract) ResultSet
    % Class that represents a set of results, for example L1O or runValidation.
    properties
        results (1,:) GMLVQ.Result
        averageRun GMLVQ.AverageRun
    end
    
    properties (Dependent)
        finalTprs
        finalFprs
        finalThresholds
        nRuns
%        averageConfusionMatrix
    end
    
    methods
        % Getter for the finalTprs
        function res=get.finalTprs(this)
            res = cell2mat(arrayfun(@(x) x.validationPerf(end).tpr, [this.results.run], 'UniformOutput', false))';
        end
        
        % Getter for the finalFprs
        function res=get.finalFprs(this)
            res = cell2mat(arrayfun(@(x) x.validationPerf(end).fpr, [this.results.run], 'UniformOutput', false))';
        end
        
        % Getter for the finalThresholds
        function res=get.finalThresholds(this)
            res = cell2mat(arrayfun(@(x) x.validationPerf(end).thresholds, [this.results.run], 'UniformOutput', false))';
        end
        
        % Getter for the nRuns
        function res=get.nRuns(this)
            res = length(this.results);
        end
    end
    
    methods (Abstract)
        % Plot the results
        plot(this);
    end
end