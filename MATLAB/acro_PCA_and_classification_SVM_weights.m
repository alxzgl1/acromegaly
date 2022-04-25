%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_PCA_and_classification_SVM_weights()

clc;

aFile = 'Lipids_N'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

bZScore = 1; % z-score

nTrainTestSplit = 0.50; % 50% (default)

pEVs = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95];
nEVs = length(pEVs);

nModels = 1000; % 5000

aModelAveraging = 'poly2'; % 'poly{1,2,3}', 'mean'

nPercentageOfSamples = 1.0; % 1.0 (default), nPercentageOfSamples *  min(patients), nPercentageOfSamples * controls
nEqualiseSamples = 1;

aAnalysisType = 'multivariate'; % KNN-imputation, norm by sum, log10, pareto-scaling

% path
aPath = 'd:\\data\\acromegaly';

% parameters
nMV_BySubjects = 0.30;
nMV_ByFeatures = 0.25; % 0.25 (default), 0.50

% load data
aExclusion = sprintf('S%dF%d', round(100 * nMV_BySubjects), round(100 * nMV_ByFeatures));
aFilename = sprintf('%s\\import\\MA\\%s\\%s_%s_NA.csv', aPath, aAnalysisType, aFile, aExclusion);
T = readcell(aFilename);
% parse data
data = cell2mat(T(3:end, 2:end)); 
subjects = T(1, 2:end);
pLabelsClass = contains(T(2, 2:end), 'Acromegaly');
names = T(3:end, 1);

% loop EV
for iEV = 1:nEVs
  nExplainedVariance = pEVs(iEV);
  % init
  Y = data; % in GMLVQ [labels x points]
  Y_labels = pLabelsClass; % [1 x labels]
  
  % svd
  [~, S, V] = svd(cov(Y')); 
  d = diag(S);
  ev = cumsum(d) / sum(d); % explained variance
  % keep K components
  th = nExplainedVariance;
  K = find(ev > th, 1);
  V = V(:, 1:K); % eig: V = V(:, i(1:K));
  % plot eigen-values
  bDebugEigenvalues = 0;
  if bDebugEigenvalues == 1
    plot(d);
  end
  % dimensionality reduction
  X = V' * Y;
  y = Y_labels; 
  
  % normalise data
  % if bZScore == 1
  %   X = (X - mean(X, 2)) ./ std(X, 0, 2);
  % end
  
  % init
  pAccuracies = zeros(nModels, 1);
  pCoeffs = zeros(K, nModels);

  % init splits
  nSplit = nTrainTestSplit;
  fprintf(1, 'split: %d/%d\n', round(100 * nSplit), round(100 * (1 - nSplit)));
  % init labels
  labels_1 = find(Y_labels == 1); 
  labels_0 = find(Y_labels == 0); 
  
  N1 = length(labels_1); labels_1 = labels_1(1:round(nPercentageOfSamples * N1));
  N0 = length(labels_0); labels_0 = labels_0(1:round(nPercentageOfSamples * N0));
  
  N1 = length(labels_1); H1 = round(nSplit * N1);
  N0 = length(labels_0); H0 = round(nSplit * N0); 
  
  if nEqualiseSamples == 1
    H = min(H1, H0); H1 = H; H0 = H;
    N = min(N1, N0); N1 = N; N0 = N;
  end
  
  % loop models with random sampling
  for iModel = 1:nModels
    fprintf(1, '%d / %d\n', iModel, nModels);
    % permute labels
    j_labels_1 = labels_1(randperm(length(labels_1)));
    j_labels_0 = labels_0(randperm(length(labels_0)));
    iTR = [j_labels_1(1:H1), j_labels_0(1:H0)];
    iTS = [j_labels_1((H1 + 1):N1), j_labels_0((H0 + 1):N0)];
    % init labels
    y_TR = y(iTR);
    y_TS = y(iTS);
    % init data
    X_TR = X(:, iTR);
    X_TS = X(:, iTS);

    % z-score
    if bZScore == 1
      X_TR = (X_TR - mean(X_TR, 2)) ./ std(X_TR, 0, 2);
      X_TS = (X_TS - mean(X_TS, 2)) ./ std(X_TS, 0, 2);
    end
    
    % SVM
    SVM_model = fitcsvm(X_TR', y_TR', 'KernelFunction', 'linear', 'KernelScale', 'auto');
    SVM_beta = SVM_model.Beta;
    % prediction
    bPredictManual = 0;
    if bPredictManual == 1 % SVM predict (manual)
      SVM_bias = SVM_model.Bias;
      SVM_scale = SVM_model.KernelParameters.Scale;
      s = (X_TS' / SVM_scale) * SVM_beta + SVM_bias;
      SVM_labels = s > 0;
    else
      SVM_labels = predict(SVM_model, X_TS'); % SVM predict (MATLAB)
    end
    % accuracy
    SVM_accuracy = mean(SVM_labels(:) == y_TS(:));
    % init
    beta = SVM_beta / vecnorm(SVM_beta);
    pAccuracies(iModel) = SVM_accuracy;
    pCoeffs(:, iModel) = beta;  
  end
 
  % save weights
  bSaveCoeffs = 1;
  if bSaveCoeffs == 1
    [w, accuracy] = support_average_models(pAccuracies, pCoeffs, aModelAveraging);
    x = w' * X;
    labels = Y_labels;
    aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\weights\\%s_%d%d_EV_%d_%s.mat', aFile, 100 * nSplit, round(100 * (1 - nSplit)), round(100 * nExplainedVariance), aModelAveraging);
    save(aFilename, 'w', 'x', 'accuracy', 'subjects', 'labels', 'K');
  end
end

end % end

%-------------------------------------------------------------------------------