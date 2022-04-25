% Function that runs L1O: it does nFeatureVectors runs with each time one feature vector removed
% from the training set.
% @param this GMLVQ.GMLVQ
% @out res GMLVQ.L1OResultSet
function res = runL1O(this)

% We continuously change some parameters, so we clone the GMLVQ object to be able to modify it
% without affecting the callee's result for future calls. This is necessary because the object is a
% handle object so it is passed by reference.
this = copy(this);
GMLVQ.Helpers.setRNG(this.params.rngseed);

nRuns = this.data.nFeatureVectors; % Each run we omit one sample
scoreval = zeros(1, this.data.nFeatureVectors); % Score of the left out example
lambda = zeros(nRuns, this.nDimensions, this.nDimensions);
protos = zeros(nRuns, this.nPrototypes, this.nDimensions);
confmat = zeros(nRuns, this.nClasses, this.nClasses);

results = GMLVQ.Result.empty(0, nRuns);

disp(['Learning curves, averages over ', num2str(nRuns), ' L1O runs']);

%ADD RJV Calculate validation confusion matrix from remaining protos?

% Do the runs
% tic
parfor krun = 1 : nRuns
    disp(['Leave one out: ', num2str(krun), ' of ', num2str(nRuns)]);
    
    trainsetind = setdiff(1:this.data.nFeatureVectors, krun); % Leave only one sample out
    fvectrain = this.data.featureVectors(trainsetind,:);
    lbltrain = this.data.labels(trainsetind);
    fvecout = this.data.featureVectors(krun,:);
    lblout = this.data.labels(krun);

    % Perform the algorithm, do a run
    trainingData = GMLVQ.DataPair(fvectrain, lbltrain);
    validationData = GMLVQ.DataPair(fvecout, lblout);
    run = GMLVQ.Run(this, trainingData, GMLVQ.DataPair.empty);
    result = run.execute();
    
    % We retrieve the score of the omitted data
    [~, crout, ~, scoreval(krun)] = result.run.classify(validationData);
    
    amat = zeros(this.nClasses, this.nClasses);
    amat(lblout,crout) = 1;
    confmat(krun, :,:) = amat; %RJV calculate confmat
    
    protos(krun, :, :) = result.run.prototypes;
    lambda(krun, :, :) = result.run.lambda;
    results(krun) = result;
end
% fprintf('[Runs] '); toc

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

% Updated mean prototypes, matrix
protos_mean = wm;
protos_std = sqrt(wm2 - wm.^2);
lambda_mean = squeeze(mean(lambda, 1));
lambda_std = sqrt(squeeze(mean(lambda.^2, 1)) - lambda_mean.^2);

% ROC data (the average already)
[tpr, fpr, auroc, thresholds] = GMLVQ.Helpers.roc(this.data.labels ~= this.params.rocClass, scoreval);

res = GMLVQ.L1OResultSet;
res.results = results;

% Construct average run
averageRun = GMLVQ.AverageRun(this, this.data);
averageRun.setMeanPerformanceTraining(results);
averageRun.prototypes = protos_mean;
averageRun.prototypesStd = protos_std;
averageRun.lambda = lambda_mean;
averageRun.lambdaStd = lambda_std;
averageRun.validationPerf(1) = GMLVQ.Performance(this.nClasses);
averageRun.validationPerf(1).tpr = tpr;
averageRun.validationPerf(1).fpr = fpr;
averageRun.validationPerf(1).auroc = auroc;
averageRun.validationPerf(1).thresholds = thresholds;
averageRun.validationPerfConfusionMatrix{1} = squeeze(sum(confmat,1));
averageRun.validationPerf(1).confusionMatrix = averageRun.validationPerfConfusionMatrix{1} / nRuns;
res.averageRun = averageRun;

end
