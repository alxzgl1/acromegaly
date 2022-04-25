classdef Performance
    % Class that stores the performance of the algorithm in terms of total error, cost function,
    % AUC(ROC) and class-wise errors.
    properties
        totalError (1,1) { mustBeNumeric, mustBeFinite } = 0
        costFunction (1,1) { mustBeNumeric, mustBeFinite } = 0
        classWise (:,1) { mustBeNumeric, mustBeFinite }
        score (:,1) { mustBeNumeric, mustBeFinite }
        auroc (1,1) { mustBeNumeric, mustBeFinite } = 0
        tpr (:,1) { mustBeNumeric, mustBeFinite }
        fpr (:,1) { mustBeNumeric, mustBeFinite }
        thresholds (:,1) { mustBeNumeric, mustBeFinite }
        confusionMatrix (:,:) { mustBeNumeric, mustBeFinite }
    end
    
    methods
        % Constructor of this class
        % @param nClasses integer
        % @out this GMLVQ.Performance
        function this = Performance(nClasses)
            % Initializes the proper size of the arrays
            this.classWise = zeros(nClasses, 1);
            this.confusionMatrix = zeros(nClasses, nClasses);
        end
    end
end

