classdef Run < matlab.mixin.Copyable
    % Class that represents one execution of the algorithm and contains all relevant data for it
    properties
        gmlvq GMLVQ.GMLVQ
        
        trainingData GMLVQ.DataPair
        trainingPerf (:,:) GMLVQ.Performance %RJV FIX Unreasonable Size constraint (1,:) -> (:,:)
        validationData GMLVQ.DataPair
        validationPerf (:,:) GMLVQ.Performance %RJV FIX Idem
        
        prototypes (:,:) { mustBeNumeric, mustBeFinite }
        lambda (:,:) { mustBeNumeric }
        meanFeatures (1,:) { mustBeNumeric }
        stdFeatures (1,:) { mustBeNumeric }
        omegaMatrix (:,:) { mustBeNumeric }
    end
    
    properties (Dependent)
        nPrototypes
        nClasses
        nDimensions
        nSteps
        
        doValidation
    end
    
    methods
        % Constructor
        % @param gmlvq GMVLQ.GMLVQ
        % @param trainingData GMLVQ.DataPair
        % @param validationData GMLVQ.DataPair
        % @out this GMLVQ.Run
        function this = Run(gmlvq, trainingData, validationData)
            % The runner modifies data here, we don't want to leak this
            % change to the user scope.
            this.gmlvq = copy(gmlvq);
            
            this.trainingData = copy(trainingData);
            this.trainingPerf = GMLVQ.Performance.empty(this.nSteps, 0);
            this.validationData = copy(validationData);
            this.validationPerf = GMLVQ.Performance.empty(this.nSteps, 0);
            
            % Calculate ZScore items
            this.meanFeatures = mean(this.trainingData.featureVectors, 1);
            this.stdFeatures = std(this.trainingData.featureVectors, 0, 1);
            
            % Apply if necessary
            if this.gmlvq.params.doztr
                this.trainingData.featureVectors = this.zScore(this.trainingData);
                if this.doValidation; this.validationData.featureVectors = this.zScore(this.validationData); end
            end
        end
        
        % Getter for the nPrototypes
        function res=get.nPrototypes(this)
            res = this.gmlvq.nPrototypes;
        end
        
        % Getter for the nClasses
        function res=get.nClasses(this)
            res = this.gmlvq.nClasses;
        end
        
        % Getter for the nDimensions
        function res=get.nDimensions(this)
            res = this.trainingData.nDimensions;
        end
        
        % Getter for the nSteps
        function res=get.nSteps(this)
            res = this.gmlvq.totalsteps;
        end
        
        % Getter for the doValidation
        function res=get.doValidation(this)
            res = ~isempty(this.validationData);
        end
        
        % Function to evaluate performance of this classifier
        % @param this GMLVQ.Run
        % @param dataPair GMLVQ.DataPair
        % @out costFunction double
        % @out crispOut 1xN double vector
        % @out margins 1xN double vector
        % @out score 1xN double vector
        % (with N = nFeatureVectors)
        function [costFunction, crispOut, margins, score] = classify(this, dataPair)
            % Check whether truth is known (labels present in testing data)
            if dataPair.nFeatureVectors ~= length(dataPair.labels)
                % Otherwise set garbage
                dataPair.labels = ones(1, dataPair.nFeatureVectors);
                dataPair.labels(ceil(dataPair.nFeatureVectors / 2), end) = 2;
            end
            
            % Optionally perform zScore transformation
            if this.gmlvq.params.doztr
                dataPair.featureVectors = this.zScore(dataPair);
            end
            
            [costFunction, crispOut, margins, score] ...
                = this.cost(dataPair, this.prototypes, sqrtm(this.lambda), 0);
        end
        
        % Function to visualise the confusion matrix
        % @param this GMLVQ.Run
        % @param idx double
        function visu_conf(this, idx) %RJV
            %   VISU_CONF visualize the confusion matrix. 
            %   requires Deep Learning Toolbox.
            %   See also: plotconfusion
            if nargin < 2
                idx = length(this.trainingPerf);
                title = 'Final';
            else
                title = sprintf('Run %d', idx);
            end
            GMLVQ.Util.confplot(this.trainingPerf(idx).confusionMatrix, title);
        end
        
        function visu_conf_ani(this) %RJV
            idx = length(this.trainingPerf);
            for i = 1:idx-1
                this.visu_conf(i);
                pause(0.1);
            end
            this.visu_conf();
        end
        
        % These functions are in other files
        plot(this);
        visu_2d(this);
        
        plotconf(this, type, stepID, suffix);
    end
    
    % These methods are only accessible by the listed classes and their
    % subclasses.
    methods (Access = { ?GMLVQ.GMLVQ, ?GMLVQ.Run })
        res = execute(this);
    end
    
    methods (Access = protected)
        % Performs a z-score transformation of fvec and returns rescaled vectors
        % @param this GMLVQ.Run
        % @param dataPair GMLVQ.DataPair
        % @out featureVectors
        function featureVectors = zScore(this, dataPair)
            featureVectors = (dataPair.featureVectors - repmat(this.meanFeatures, dataPair.nFeatureVectors, 1)) ...
                    ./ repmat(this.stdFeatures, dataPair.nFeatureVectors, 1);
        end
        
        % Calculates the performance of the given data with the given prototypes and omega matrix,
        % using cost and roc functions.
        % @param this GMLVQ.Run
        % @param dataPair GMLVQ.DataPair
        % @param prototypes
        % @param omegaMatrix
        % @param customMu double optional
        % @out performance GMLVQ.Performance
        function performance = calculatePerformance(this, dataPair, prototypes, omegaMatrix, customMu)
            performance = GMLVQ.Performance(this.nClasses);
            
            % Calculate costs
            if nargin < 5
                [performance.costFunction, crispOut, margins, performance.score] ...
                    = this.cost(dataPair, prototypes, omegaMatrix);
            else
                [performance.costFunction, crispOut, margins, performance.score] ...
                    = this.cost(dataPair, prototypes, omegaMatrix, customMu);
            end
            performance.totalError = sum(margins > 0) / dataPair.nFeatureVectors;
            
            % Calculate ROC
            [performance.tpr, performance.fpr, performance.auroc, performance.thresholds]...
                = GMLVQ.Helpers.roc(dataPair.labels ~= this.gmlvq.params.rocClass, performance.score);
            
            % Classwise errors
            classWise = zeros(1, this.nClasses);
            for i = 1:this.nClasses
                classWise(i) = sum(margins(dataPair.labels == i) > 0) / sum(dataPair.labels == i);
            end
            performance.classWise = classWise;
            
            % Confusion matrix
            for i = 1:dataPair.nFeatureVectors
                performance.confusionMatrix(dataPair.labels(i), crispOut(i))...
                    = performance.confusionMatrix(dataPair.labels(i), crispOut(i)) + 1;
            end
        end
        
        % These functions are in other files
        [prototypes, omegaMatrix] = setInitialPrototypes(this);
        [costFunction, crispOut, margins, score] = cost(this, dataPair, prototypes, omegaMatrix, customMu);
        [prototypes, omegaMatrix] = doBatchStep(this, prototypes, omegaMatrix);
    end    
end