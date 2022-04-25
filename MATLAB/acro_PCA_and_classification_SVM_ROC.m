%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_PCA_and_classification_SVM_ROC()

clc;

bShuffled = 0; % test

aFile = 'HILIC_P'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

bZScore = 1; % z-score

nEV = 0.75;

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
th = nEV;
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

% init labels
labels_1 = find(Y_labels == 1); 
labels_0 = find(Y_labels == 0); 

% N1 = length(labels_1); labels_1 = labels_1(1:N1);
% N0 = length(labels_0); labels_0 = labels_0(1:N0);

N1 = length(labels_1); % H1 = N1;
N0 = length(labels_0); % H0 = N0; 

if nEqualiseSamples == 1
  N = min(N1, N0); N1 = N; N0 = N;
end

% permute labels
% j_labels_1 = labels_1(randperm(length(labels_1)));
% j_labels_0 = labels_0(randperm(length(labels_0)));
% iTR = [j_labels_1(1:N1), j_labels_0(1:N0)];

iTR = [labels_1(1:N1), labels_0(1:N0)];

% init labels
y_TR = y(iTR);
% init data
X_TR = X(:, iTR);

% shuffle labels
if bShuffled == 1
  y_TR = y_TR(randperm(length(y_TR)));
end

% z-score
if bZScore == 1
  X_TR = (X_TR - mean(X_TR, 2)) ./ std(X_TR, 0, 2);
end

% SVM
mdl = fitcsvm(X_TR', y_TR', 'KernelFunction', 'linear', 'KernelScale', 'auto');

% compute the posterior probabilities (scores)
mdl = fitPosterior(mdl);
[~, scores] = resubPredict(mdl);

% compute the standard ROC curve using the scores from the SVM model.
% [svm_X, svm_Y, svm_T, svm_AUC] = perfcurve(y_TR', scores(:, mdl.ClassNames), 1);
[svm_X, svm_Y, svm_T, svm_AUC] = perfcurve(y_TR', scores(:, 2), 1);

% display the area under the curve
svm_AUC

% plot the ROC curve
plot(svm_X, svm_Y);
xlabel('False positive rate');
ylabel('True positive rate');
  
% accuracy
% fprintf(1, '%s: accuracy = %1.2f\n', aEnsemble, 100 * mean(y_ens == Y_test));

end % end

%-------------------------------------------------------------------------------