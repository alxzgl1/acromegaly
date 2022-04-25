%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_PCA_and_classification_SVM_ROC()

clc;

aFile = 'HILIC_P'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

bZScore = 1; % z-score

nTrainTestSplit = 1.0; 

nExplainedVariance = 0.75;

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

% loop shuffles
for bShuffled = 0:1
  % init splits
  nSplit = nTrainTestSplit;
  fprintf(1, 'split: %d/%d | shuffled (%d)\n', round(100 * nSplit), round(100 * (1 - nSplit)), bShuffled);
  % init labels
  labels_1 = find(Y_labels == 1); 
  labels_0 = find(Y_labels == 0); 

  % labels
  labels_1 = labels_1(randperm(length(labels_1)));
  labels_0 = labels_0(randperm(length(labels_0)));

  N1 = length(labels_1); H1 = round(nSplit * N1);
  N0 = length(labels_0); H0 = round(nSplit * N0); 

  if nEqualiseSamples == 1
    H = min(H1, H0); H1 = H; H0 = H;
    N = min(N1, N0); N1 = N; N0 = N;
  end

  iTR = [labels_1(1:H1), labels_0(1:H0)];
  iTS = [labels_1((H1 + 1):N1), labels_0((H0 + 1):N0)];

  % init labels
  y_TR = y(iTR);
  y_TS = y(iTS);

  % init data
  X_TR = X(:, iTR);
  X_TS = X(:, iTS);

  % shuffle labels
  if bShuffled == 1
    y_TR = y_TR(randperm(length(y_TR)));
    y_TS = y_TS(randperm(length(y_TS)));
  end

  % z-score
  if bZScore == 1
    X_TR = (X_TR - mean(X_TR, 2)) ./ std(X_TR, 0, 2);
    X_TS = (X_TS - mean(X_TS, 2)) ./ std(X_TS, 0, 2);
  end

  % init ensemble learning
  X_train = X_TR';
  Y_train = y_TR';
  X_test = X_TS';
  Y_test = y_TS';

  % i_train = randperm(length(Y_train));
  % i_test = randperm(length(Y_test));

  % X_train = X_train(i_train, :);
  % Y_train = Y_train(i_train, :);

  % X_test = X_test(i_test, :);
  % Y_test = Y_test(i_test, :);

  % SVM
  mdl = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'KernelScale', 'auto');
  % compute the posterior probabilities (scores)
  mdl = fitPosterior(mdl);
  [~, scores] = resubPredict(mdl);
  % compute the standard ROC curve using the scores from the SVM model.
  % [svm_X, svm_Y, svm_T, svm_AUC] = perfcurve(Y_train, scores(:, mdl.ClassNames), 1);
  [svm_X, svm_Y, svm_T, svm_AUC] = perfcurve(Y_train, scores(:, 2), 1);

  % display the area under the curve
  svm_AUC
  
  % plot the ROC curve
  plot(svm_X, svm_Y);
  xlabel('False positive rate');
  ylabel('True positive rate');
    
  % accuracy
  % fprintf(1, '%s: accuracy = %1.2f\n', aEnsemble, 100 * mean(y_ens == Y_test));
end

end % end

%-------------------------------------------------------------------------------