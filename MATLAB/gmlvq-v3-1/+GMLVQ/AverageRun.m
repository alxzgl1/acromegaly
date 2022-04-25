classdef AverageRun < GMLVQ.Run
    % Class that represents the "average" of a number of Runs, for example in a L1O-scenario.
    
    properties
        trainingPerfStd (:,:) GMLVQ.Performance %FIX RJV
        validationPerfStd (:,:) GMLVQ.Performance %FIX RJV
        
        prototypesStd (:,:) { mustBeNumeric, mustBeFinite }
        lambdaStd (:,:) { mustBeNumeric }

        trainingPerfConfusionMatrix cell % ADD RJV Summed confusionMatrices        
        validationPerfConfusionMatrix cell      
    end
    
    methods
        % Constructor of this class
        % @param gmlvq GMLVQ.GMLVQ
        % @param trainingData GMLVQ.DataPair
        % @out GMLVQ.AverageRun
        function this = AverageRun(gmlvq, trainingData)
            this@GMLVQ.Run(gmlvq, trainingData, GMLVQ.DataPair.empty);
            
            this.trainingPerfStd = GMLVQ.Performance.empty(this.nSteps, 0);
            this.validationPerfStd = GMLVQ.Performance.empty(this.nSteps, 0);
        end
        
        % Function that averages the validation performance parameters of the given result set
        % @param this GMLVQ.AverageRun
        % @param results GMLVQ.Result[]
        function setMeanPerformanceTraining(this, results)
            arguments
                this GMLVQ.AverageRun
                results (:,1) GMLVQ.Result
            end
            this.setMeanPerformance(results, 'trainingPerf');
        end
        
        % Function that averages the training performance parameters of the given result set
        % @param this GMLVQ.AverageRun
        % @param results GMLVQ.Result[]
        function setMeanPerformanceValidation(this, results)
            arguments
                this GMLVQ.AverageRun
                results (:,1) GMLVQ.Result
            end
            this.setMeanPerformance(results, 'validationPerf');
        end

        % Function to visualise the confusion matrix
        % @param this GMLVQ.Run
        % @param idx double
        function f = visu_conf(this, idx) %RJV
            if nargin < 3
                idx = -1;
            end
            
            if isempty(this.validationPerfConfusionMatrix) % RJV Workaround for L1O
                f = figure;
                this.visu_conf_select('t', idx);
                scrsz = get(0,'ScreenSize');
                f.Position(4) = round(scrsz(4)/2.5); %make height = width / 2
                f.Position(3) = f.Position(4);
            else                
                f = figure('Name', 'Confusion Matrices');
                scrsz = get(0,'ScreenSize');
                f.Position(3) = round(scrsz(3) / 2.8);
                f.Position(4) = round(f.Position(3) / 2.1); %make height = width / 2
                ax1 = subplot(1,2,1,'parent',f);
                ax2 = subplot(1,2,2,'parent',f);
                
                figure; % plotconf will use gcf so prevent overwriting
                h = this.visu_conf_select('t', idx);
                
                GMLVQ.Util.cmcToSubplot(h, f, ax1);
                
                figure; % plotconf will use gcf so prevent overwriting
                h = this.visu_conf_select('f', idx);
                GMLVQ.Util.cmcToSubplot(h, f, ax2);
            end
        end
        
        % Function to visualise the confusion matrix
        % @param this GMLVQ.Run
        % @param idx double
        function h = visu_conf_select(this, type, idx) %RJV
            %   VISU_CONF visualize the confusion matrix. 
            %   type 't' or 'f'
            %   requires Deep Learning Toolbox.
            %   See also: plotconfusion
            if nargin < 3 || idx == -1
                title = 'Final';
            else
                title = sprintf('Run %d', idx);
            end
            if nargin < 2
                type = 't';
            end
            if type == 't'
                type = 'trainingPerfConfusionMatrix';
                title = sprintf('%s (Training) \\Sigma', title);                
            else
                type = 'validationPerfConfusionMatrix';
                title = sprintf('%s (Validation) \\Sigma', title);
            end
            if nargin < 3 || idx == -1
                idx = length(this.(type));
            end
            if ~isempty(this.(type))
                h = GMLVQ.Util.confplot(this.(type){idx}, title);
            else
                warning off backtrace
                warning(['No ConfusionMatrix for [', type, '], results unpredictable']);
                warning on backtrace 
                h = figure;
                axes(h);
            end
        end        
        
    end
    
    methods (Access = protected)
        % Function that averages the performance parameters of the given result set
        % @param this GMLVQ.AverageRun
        % @param results GMLVQ.Result[]
        % @param type string
        function setMeanPerformance(this, results, type)
            arguments
                this GMLVQ.AverageRun
                results (:,1) GMLVQ.Result
                type (1,:) char
            end
            % Each result contains a list of Performance objects that records the performance at a
            % specific step. We want to average the various metrics across the results, but on a
            % step-by-step basis. So the performance of run 1, step 1 is averaged with that of run
            % 2, step 1; run 3, step 1; etc.
            
            stdType = [type 'Std'];
            nResults = length(results);
            this.(type)(1:this.nSteps + 1) = GMLVQ.Performance(this.nClasses);
            this.(stdType)(1:this.nSteps + 1) = GMLVQ.Performance(this.nClasses);
            
            % This construct puts all auroc-vectors side by side. It makes a list of all runs, and
            % for each run it takes all aurocs in a column vector (the transpose)
            % Using cell magic we can then easily assign them to the separate step-performances
            aurocs = cell2mat(arrayfun(@(x) [x.(type).auroc]', [results.run], 'UniformOutput', false));
            aurocMean = num2cell(mean(aurocs, 2));
            aurocStd = num2cell(std(aurocs, 1, 2)); % We normalize by N, not N-1
            [this.(type).auroc] = aurocMean{:};
            [this.(stdType).auroc] = aurocStd{:};
            
            % This construct puts all classWise-matrices side by side in the third dimension
            classWiseCells = arrayfun(@(x) [x.(type).classWise]', [results.run], 'UniformOutput', false);
            classWiseSize = [(this.nSteps + 1) this.nClasses nResults]; % Take the size of the separate cw arrays
            classWise = reshape(cell2mat(classWiseCells), classWiseSize); % Put the matrices in another array
            classWiseMean = num2cell(squeeze(mean(classWise, 3))', 1);
            classWiseStd = num2cell(squeeze(std(classWise, 1, 3))', 1);
            [this.(type).classWise] = classWiseMean{:};
            [this.(stdType).classWise] = classWiseStd{:};
            
            % This construct puts all costFunction-vectors side by side.
            costFunctions = cell2mat(arrayfun(@(x) [x.(type).costFunction]', [results.run], 'UniformOutput', false));
            costFunctionMean = num2cell(mean(costFunctions, 2));
            costFunctionStd = num2cell(std(costFunctions, 1, 2));
            [this.(type).costFunction] = costFunctionMean{:};
            [this.(stdType).costFunction] = costFunctionStd{:};
            
            % This construct puts all totalError-vectors side by side.
            totalErrors = cell2mat(arrayfun(@(x) [x.(type).totalError]', [results.run], 'UniformOutput', false));
            totalErrorMean = num2cell(mean(totalErrors, 2));
            totalErrorStd = num2cell(std(totalErrors, 1, 2));
            [this.(type).totalError] = totalErrorMean{:};
            [this.(stdType).totalError] = totalErrorStd{:};
            
            % This construct puts all tpr-matrices side by side in the third dimension
            tprCells = arrayfun(@(x) [x.(type).tpr]', [results.run], 'UniformOutput', false);
            tprSize = [(this.nSteps + 1) size(tprCells{1}, 2) nResults];
            tpr = reshape(cell2mat(tprCells), tprSize);
            tprMean = num2cell(squeeze(mean(tpr, 3))', 1);
            tprStd = num2cell(squeeze(std(tpr, 1, 3))', 1);
            [this.(type).tpr] = tprMean{:};
            [this.(stdType).tpr] = tprStd{:};
            
            % This construct puts all fpr-matrices side by side in the third dimension
            fprCells = arrayfun(@(x) [x.(type).fpr]', [results.run], 'UniformOutput', false);
            fprSize = [(this.nSteps + 1) size(fprCells{1}, 2) nResults];
            fpr = reshape(cell2mat(fprCells), fprSize);
            fprMean = num2cell(squeeze(mean(fpr, 3))', 1);
            fprStd = num2cell(squeeze(std(fpr, 1, 3))', 1);
            [this.(type).fpr] = fprMean{:};
            [this.(stdType).fpr] = fprStd{:};

            confType = [type 'ConfusionMatrix'];
            confCells = zeros(this.nClasses, this.nClasses, nResults);
            this.(confType) = cell(1,this.nSteps+1);
%             tic
            for rno = 1:this.nSteps+1
                for res = 1:nResults
                    confCells(:,:,res) = results(res).run.(type)(rno).confusionMatrix;
                end
                this.(type)(rno).confusionMatrix = mean(confCells,3); 
                this.(confType){rno} = sum(confCells,3); 
                this.(stdType)(rno).confusionMatrix = std(confCells,1,3); 
                % TODO: std sometimes gives strange results, but they seem to check out.                
            end
%             fprintf('[ConfMat] ');
%             toc
%             this.(type).confusionMatrix = mean(confCells,4);
%             this.(stdType).confusionMatrix = std(confCells,1,4);

            
            % This construct puts all confusion-matrices side by side in the third dimension
%             confCells = arrayfun(@(x) [x.(type).confusionMatrix]', [results.run], 'UniformOutput', false);
%             confSize = [(this.nSteps + 1) size(confCells{1}, 2) nResults];
%             conf = reshape(cell2mat(confCells), confSize);
%             confMean = num2cell(squeeze(mean(conf, 3))', 1);
%             confStd = num2cell(squeeze(std(conf, 1, 3))', 1);
%             [this.(type).confusionMatrix] = confMean{:};
%             [this.(stdType).confusionMatrix] = confStd{:};
            
%             ConfusionCat = results(1).run.(type).confusionMatrix;            
%             for i = 2:nResults %RJV: Could be slightly more efficient.
%               ConfusionCat = cat(3, ConfusionCat, results(i).run.(type).confusionMatrix);
%             end
%             
%             conMean = mean(ConfusionCat,3);
%             [this.(type).confusionMatrix] = conMean{:};
%             conStd = std(ConfusionCat,1,3);
%             [this.(stdType).confusionMatrix] = conStd{:};            
            
            %RJV: Added missing Confusion Matrix Construction %TODO!!           
            % Reference code from 2.4: run_validation.m:138
            
%             results(1).run.trainingPerf.confusionMatrix
            
            % this is the confusion matrix for each run, we have that
%             for ios = 1:length(lblout)  % loop through all classes
%                 % computation of actual confusion matrix
%                 confact(lblout(ios),crout(ios))=confact(lblout(ios),crout(ios))+1;
%             end
            % !! plotconfusion
            % we have
%             results(1).run.trainingPerf.confusionMatrix
%             results(1).run.validationPerf.confusionMatrix
            % so, in this function, that will be:
% %             results(1).run.(type).confusionMatrix
            
            %now, for each run, plug it in here: %maybe parfor?
            % calculate averaged confusion matrix (percentages) over runs
%             for icr = 1:size(confmat,1)
%                 confmat(icr,:)=confmat(icr,:)+confact(icr,:)/sum(confact(icr,:))*100/nruns;
%             end
            
        end
    end
end