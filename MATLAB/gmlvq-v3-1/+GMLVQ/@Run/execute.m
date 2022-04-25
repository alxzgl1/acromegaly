% Function that performs the algorithm for the current Run instance
% @param this GMLVQ.Run
% @out res GMLVQ.Result
function res = execute(this)

stepsizeMatrix = zeros(this.nSteps + 1, 1);
stepsizePrototypes = zeros(this.nSteps + 1, 1);

[prototypes, omegaMatrix] = this.setInitialPrototypes();
% Copies of prototypes and omegas are stored for the adaptive step size procedure later on
prototypeCopies = zeros(this.gmlvq.params.ncop, size(prototypes, 1), size(prototypes, 2));
omegaCopies = zeros(this.gmlvq.params.ncop, size(omegaMatrix, 1), size(omegaMatrix, 2));

% We calculate the initial values of the performance for the learning curves
this.trainingPerf(1) = this.calculatePerformance(this.trainingData, prototypes, omegaMatrix);
if this.doValidation; this.validationPerf(1) = this.calculatePerformance(this.validationData, prototypes, omegaMatrix); end
stepsizeMatrix(1) = this.gmlvq.params.etam;
stepsizePrototypes(1) = this.gmlvq.params.etap;

% Perform the first (ncop) steps of the gradient descent
for step = 1 : this.gmlvq.params.ncop
    % Actual work being done
    [prototypes, omegaMatrix] = this.doBatchStep(prototypes, omegaMatrix);
    prototypeCopies(step, :, :) = prototypes;
    omegaCopies(step, :, :) = omegaMatrix;
    stepsizeMatrix(step + 1) = this.gmlvq.params.etam;
    stepsizePrototypes(step + 1) = this.gmlvq.params.etap;
    
    % Determine performance
    if this.doValidation % validation has mu=0 for the performance TODO HWY
        omegaMatrix = omegaMatrix / norm(omegaMatrix, 'fro'); % TODO WHY?!
        
        % Performance for the training data
        this.trainingPerf(step + 1) = this.calculatePerformance(this.trainingData, prototypes, omegaMatrix, 0);
    
        % Performance for the validation data
        this.validationPerf(step + 1) = this.calculatePerformance(this.validationData, prototypes, omegaMatrix, 0);
    else
        % Performance for the training data only
        this.trainingPerf(step + 1) = this.calculatePerformance(this.trainingData, prototypes, omegaMatrix, 0);
    end
end

% Compute further steps with the dynamic step sizes
for step = (this.gmlvq.params.ncop + 1):this.gmlvq.totalsteps
    % Calculate mean positions over latest steps
    protmean = squeeze(mean(prototypeCopies, 1));
    ommean = squeeze(mean(omegaCopies, 1));
    ommean = ommean / norm(ommean, 'fro'); % normalization does not change results
    
    % Compute cost mean prototypes or mean omega
    costmp = this.cost(this.trainingData, protmean, omegaMatrix);
    costmm = this.cost(this.trainingData, prototypes, ommean);
    
    % Cache old positions for Papari procedure
    omegaPrev = omegaMatrix;
    protPrev = prototypes;
    
    % Perform next step
    [prototypes, omegaMatrix] = this.doBatchStep(prototypes, omegaMatrix);
    
    % By default, step sizes are increased in every step, to enforce oscillatory behaviour
    this.gmlvq.params.etam = this.gmlvq.params.etam * this.gmlvq.params.incfac;
    this.gmlvq.params.etap = this.gmlvq.params.etap * this.gmlvq.params.incfac;
    
    % Cost function values to compare to old ones for the Papari procedure
    costfp = this.cost(this.trainingData, prototypes, omegaPrev, 0);
    costfm = this.cost(this.trainingData, protPrev, omegaMatrix);
    
    % Heuristic extension of the Papari procedure
    % Treats matrix and prototype step sizes separately
    if costmp <= costfp
        this.gmlvq.params.etap = this.gmlvq.params.etap / this.gmlvq.params.decfac;
        prototypes = protmean;
    end
    if costmm <= costfm
        this.gmlvq.params.etam = this.gmlvq.params.etam / this.gmlvq.params.decfac;
        omegaMatrix = ommean;
    end
    
    % Update the copies of the latest steps
    prototypeCopies = circshift(prototypeCopies, -1, 1);
    omegaCopies = circshift(omegaCopies, -1, 1);
    stepsizeMatrix(step + 1) = this.gmlvq.params.etam;
    stepsizePrototypes(step + 1) = this.gmlvq.params.etap;
    
    % TODO make end
    prototypeCopies(this.gmlvq.params.ncop, :, :) = prototypes;
    omegaCopies(this.gmlvq.params.ncop, :, :) = omegaMatrix;
    
    % Calculate performance
    % Performance for the training data
    this.trainingPerf(step + 1) = this.calculatePerformance(this.trainingData, prototypes, omegaMatrix, 0);

    % Performance for the validation data
    if this.doValidation
        this.validationPerf(step + 1) = this.calculatePerformance(this.validationData, prototypes, omegaMatrix, 0);
    end
end

this.lambda = omegaMatrix'*omegaMatrix;
this.omegaMatrix = omegaMatrix;
this.prototypes = prototypes;

% Construct result
res = GMLVQ.Result;
res.run = this;
res.stepsizeMatrix = stepsizeMatrix;
res.stepsizePrototypes = stepsizePrototypes;

end