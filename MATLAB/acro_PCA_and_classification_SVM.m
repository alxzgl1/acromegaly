%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_PCA_and_classification_SVM()

clc;

aFile = 'HILIC_N'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

bZScore = 1; % z-score

nTrainTestSplit = 0.50; % 50% (default)

nExplainedVariance = 0.75;

nModels = 1000; % 5000

aModelAveraging = 'poly2'; % 'poly{1,2,3}', 'mean'

bSaveImportantMetabolites = 0;

bKeepImportantMetabolites = 0;
bNumberOfImportantMetabolites = 0;

bRemoveImportantMetabolites = 0;

bKMeans = 0;

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

% init
Y = data; % in GMLVQ [labels x points]
Y_labels = pLabelsClass; % [1 x labels]

% keep metabolites
if bKeepImportantMetabolites == 1
  nSplit = nTrainTestSplit;
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\metabolites\\%s_%d%d.mat', aFile, 100 * nSplit, round(100 * (1 - nSplit)));
  load(aFilename, 'pMetaboliteIndices');
  nMetaboliteIndices = bNumberOfImportantMetabolites;
  pMetaboliteIndices = pMetaboliteIndices((end - nMetaboliteIndices + 1):end);
  Y = Y(pMetaboliteIndices, :);
end

if bRemoveImportantMetabolites > 0
  nSplit = nTrainTestSplit;
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\metabolites\\%s_%d%d.mat', aFile, 100 * nSplit, round(100 * (1 - nSplit)));
  load(aFilename, 'pMetaboliteIndices');
  nMetaboliteIndices = bRemoveImportantMetabolites;
  % keep less important
  % pMetaboliteIndices = pMetaboliteIndices(1:(end - nMetaboliteIndices));
  % keep more important
  pMetaboliteIndices = pMetaboliteIndices((end - nMetaboliteIndices):end);
  Y = Y(pMetaboliteIndices, :);
end

% k-means
if bKMeans == 1
  [~, C] = kmeans(Y, 50);
  Y = C;
end

% models
if bKeepImportantMetabolites == 1
  nModels = 1;
end

% debug 
% bDebugSignals = 0;
% if bDebugSignals == 1
%   Y = zeros(3000, 62);
%   t = linspace(0, 0.1, 62);
%   for i = 1:100
%     y = sin(2 * pi * t * i) + 0.01 * randn(1, length(t));
%     Y(i, :) = y;
%   end
%   Y(101:end, :) = 0.25 * randn(3000 - 100, length(t));
%   lot(t, Y);
% end
% debug^

% eig
% [V, D] = eig(cov(Y'));
% [d, i] = sort(diag(D), 'descend');
% x = V' * Y;

% PCA
if bKeepImportantMetabolites == 0 && bKMeans == 0
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
else
  nMetabolites = size(Y, 1);
  X = Y;
  y = Y_labels; 
  V = eye(nMetabolites);
  S = eye(nMetabolites);
  K = nMetabolites;
end

% normalise data
% if bZScore == 1
%   X = (X - mean(X, 2)) ./ std(X, 0, 2);
% end

% init
nMetabolites = size(V, 1);
pAccuracies = zeros(nModels, 2);
pCoeffs = zeros(K, nModels, 2);
pProjections = zeros(nMetabolites, nModels, 2);

pMetaboliteIndices = zeros(nMetabolites, 1); % indices

nSubjects = size(X, 2);
pXY = zeros(nSubjects, 2);

% loop shuffles
for bShuffled = 0:1
  % init splits
  nSplit = nTrainTestSplit;
  fprintf(1, 'split: %d/%d | shuffled (%d)\n', round(100 * nSplit), round(100 * (1 - nSplit)), bShuffled);
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
    % shuffle
    if bShuffled == 1
      y_TR = y_TR(randperm(length(y_TR)));
      y_TS = y_TS(randperm(length(y_TS)));
    end

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
    pAccuracies(iModel, bShuffled + 1) = SVM_accuracy;
    pCoeffs(:, iModel, bShuffled + 1) = beta; 
    pProjections(:, iModel, bShuffled + 1) = V * S(1:K, 1:K) * beta; 
  end
  % histogram
  pBins = 0.0:0.1:1;
  subplot(2, 3, 1);
  h = histc(pAccuracies(:, bShuffled + 1), pBins); hold on;
  plot(pBins, h / nModels, 'o-'); xlim([0.25, 0.95]); box off; 
  xlabel('accuracy'); ylabel('counts / sum');
  title(sprintf('hist | train/test: %d/%d | %s', 100 * nSplit, round(100 * (1 - nSplit)), aFile), 'FontWeight', 'normal');

  % coeffs
  subplot(2, 3, 2);
  w = support_average_models(pAccuracies(:, bShuffled + 1), pCoeffs(:, :, bShuffled + 1), aModelAveraging);
  plot(w, 'o-'); hold on; 
  xlim([0.5, K + 0.5]); box off; 
  xlabel('PC index'); ylabel('weights');
  title(sprintf('w | K=%d | train/test: %d/%d | %s', K, 100 * nSplit, round(100 * (1 - nSplit)), aModelAveraging), 'FontWeight', 'normal');

  % separation
  w = support_average_models(pAccuracies(:, bShuffled + 1), pCoeffs(:, :, bShuffled + 1), aModelAveraging);
  x = w' * X;
  yMin = min(x);
  yMax = max(x);
  pXY(:, bShuffled + 1) = x;
  subplot(2, 3, 3);
  if bShuffled == 0
    fill([1, 1, N1, N1], [yMin, yMax, yMax, yMin], [1.0, 0.9, 0.9], 'LineStyle', 'none'); hold on;
    line([1, nSubjects], [0, 0], 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2); hold on;
  end
  plot(x, 'o-'); hold on; 
  if bShuffled == 0
    box off; xlim([1, size(X, 2)]); ylim([yMin, yMax]);
    xlabel('patients / controls'); ylabel('w^T * X');
  end
  if bShuffled == 1
    XY0 = pXY(:, 1) > 0; 
    XY1 = pXY(:, 2) > 0; 
    ACC_XY0 = round(100 * mean(XY0(:) == Y_labels(:)));
    ACC_XY1 = round(100 * mean(XY1(:) == Y_labels(:)));
    title(sprintf('w^T * X | acc: %d%%/%d%%', ACC_XY0, ACC_XY1), 'FontWeight', 'normal');
  end
  
  % projections
  subplot(2, 3, [4, 6]);
  w = support_average_models(pAccuracies(:, bShuffled + 1), pCoeffs(:, :, bShuffled + 1), aModelAveraging);
  x = abs(V * S(1:K, 1:K) * w); % abs(V * w)
  if bShuffled == 0
    [~, ix] = sort(x);
    pMetaboliteIndices(:, 1) = ix;
  end
  x = x(ix);

  % title
  plot(x); hold on; box off; xlim([0, length(V)]);
  title(sprintf('V * w | train/test: %d/%d | PCA: %d', 100 * nSplit, round(100 * (1 - nSplit)), 1 - bKeepImportantMetabolites), 'FontWeight', 'normal');
  xlabel('metabolites'); ylabel('V * w');
end

fprintf(1, 'K = %d\n', K);

% save strongest metabolites
if bSaveImportantMetabolites == 1
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\metabolites\\%s_%d%d_EV_%d.mat', aFile, 100 * nSplit, round(100 * (1 - nSplit)), round(100 * nExplainedVariance));
  save(aFilename, 'pMetaboliteIndices', 'names');
end

% save weights
bSaveCoeffs = 0;
if bSaveCoeffs == 1
  w = squeeze(mean(pCoeffs(:, :, 1), 2));
  x = w' * X;
  labels = Y_labels;
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\weights\\%s_%d%d.mat', aFile, 100 * nSplit, round(100 * (1 - nSplit)));
  save(aFilename, 'w', 'x', 'subjects', 'labels');
end

end % end

%-------------------------------------------------------------------------------