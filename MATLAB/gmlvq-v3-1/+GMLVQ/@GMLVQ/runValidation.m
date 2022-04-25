% This function runs a validation set on the given data. This consists of [nRuns] runs in which each
% time [percentage]% of the data is left out per class and used to validate/test the obtained 
% classifier. In the end, the mean classifier is taken as result.
% @param this GMLVQ.GMLVQ
% @param nRuns integer
% @param percentage double
% @param stratify boolean
% @out res GMLVQ.ValidationResultSet
function res = runValidation(this, nRuns, percentage)

% We continuously change some parameters, so we clone the GMLVQ object to be able to modify it
% without affecting the callee's result for future calls. This is necessary because the object is a
% handle object so it is passed by reference.
this = copy(this);

% TODO make use of the RandStreams, that might be neater
GMLVQ.Helpers.setRNG(this.params.rngseed);

% Set defaults if necessary
if nargin < 3 || isempty(percentage); percentage = 10; end
if nargin < 2 || isempty(nRuns); nRuns = 5; end

disp(['Learning curves, averages over ', num2str(nRuns), ' validation runs']);
disp(['with ', num2str(percentage), '% of examples per class left out for testing']);

% Initialize all observed quantities
lambda = zeros(nRuns, this.nDimensions, this.nDimensions);
protos = zeros(nRuns, this.nPrototypes, this.nDimensions);

results = GMLVQ.Result.empty(0, nRuns);

% Do the validation runs
parfor (krun = 1 : nRuns) % set , 0 for serial work or , inf for default
    disp(['Validation run ', num2str(krun), ' of ', num2str(nRuns)]);
    
    % We create a random partition of the data into a testing and training set. The stratify option
    % indicates whether all classes should be equally represented (relative) or not.
    partition = cvpartition(this.data.labels, 'HoldOut', percentage / 100, 'Stratify', true); %!ok<PFBNS>
    
    trainingLabels = this.data.labels(partition.training(1));
    trainingFV = this.data.featureVectors(partition.training(1),:);
    validationLabels = this.data.labels(partition.test(1));
    validationFV = this.data.featureVectors(partition.test(1),:);
    
    % Perform the algorithm: do a run
    trainingData = GMLVQ.DataPair(trainingFV, trainingLabels);
    validationData = GMLVQ.DataPair(validationFV, validationLabels);
    run = GMLVQ.Run(this, trainingData, validationData);
    result = run.execute();
    
    protos(krun, :, :) = result.run.prototypes;
    lambda(krun, :, :) = result.run.lambda;
    
    results(krun) = result;
end

% Sort and average prototypes properly
wm = squeeze(protos(1,:,:)) / nRuns;
wm2 = nRuns * wm.^2;

for ipt = 2:nRuns
    protlist = 1:this.nPrototypes;
    for i = 1:this.nPrototypes
        dij = NaN(1, this.nPrototypes);
        for j = protlist
            if this.plbl(j) == this.plbl(i)
                dij(j) = norm((wm(i,:)' - squeeze(protos(ipt, j, :))), 2);
            end
        end
        
        [~, jmin] = min(dij);
        
        % Update mean and mean square
        wm(i,:) = wm(i,:) + (squeeze(protos(ipt, jmin, :))') / nRuns;
        wm2(i,:) = wm2(i,:) + (squeeze(protos(ipt, jmin, :))').^2 / nRuns;
        protlist = setdiff(protlist, jmin);
    end
end

% Updated mean prototypes, matrix, cost function, training errors, auc and class-wise errors
protos_mean = wm;
protos_std = sqrt(wm2 - wm.^2);
lambda_mean = squeeze(mean(lambda, 1));
lambda_std = sqrt(squeeze(mean(lambda.^2, 1)) - lambda_mean.^2);

res = GMLVQ.ValidationResultSet;
res.results = results;

% Construct result
averageRun = GMLVQ.AverageRun(this, this.data);
averageRun.setMeanPerformanceTraining(results);
averageRun.setMeanPerformanceValidation(results);
averageRun.prototypes = protos_mean;
averageRun.prototypesStd = protos_std;
averageRun.lambda = lambda_mean;
averageRun.lambdaStd = lambda_std;

res.averageRun = averageRun;

end
