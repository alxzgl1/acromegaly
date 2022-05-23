%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_PCA_and_classification_GMLVQ_3x()

clc;

bShuffle = 0; % test

% add toolbox
addpath('gmlvq-v3-1');

aFile = 'HILIC_N'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

bZScore = 1;

nEqualiseSamples = 0;

aAnalysisType = 'multivariate'; % KNN-imputation, norm by sum, log10, pareto-scaling

nEV = 0.75; % explained variance

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

% svd
[~, S, V] = svd(cov(Y')); 
d = diag(S);

% explained variance
ev = cumsum(d) / sum(d);

% keep K components
th = nEV;
K = find(ev > th, 1);
V = V(:, 1:K); % eig: V = V(:, i(1:K));

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

% shuffle
if bShuffle == 1
  y_TR = y_TR(randperm(length(y_TR)));
end

% z-score
if bZScore == 1
  X_TR = (X_TR - mean(X_TR, 2)) ./ std(X_TR, 0, 2);
end

% validation
GMLVQ_model = GMLVQ.GMLVQ(X_TR', (y_TR + 1)', GMLVQ.Parameters(), 30); % 30 == totalsteps
result = GMLVQ_model.runValidation(10, 25); % 10 validation runs, 25% of examples per class
lambda = result.averageRun.lambda;
prototypes = result.averageRun.prototypes;

% plot

subplot(2, 2, 1); 
plot(result.averageRun.validationPerf(end).fpr', result.averageRun.validationPerf(end).tpr', ...
  'LineWidth', 2); box off;
hold on;
% plot(mean(result.finalFprs(:, result.finalThresholds(end, :) == 0.5)), ...
%   mean(result.finalTprs(:, result.finalThresholds(end, :) == 0.5)), 'ko', ... 
%   'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot([0 1], [0 1], 'k:'); 
title(sprintf('AUC = %1.2f', -trapz(result.averageRun.validationPerf(end).fpr, ...
  result.averageRun.validationPerf(end).tpr)), 'FontWeight', 'normal');
xlabel('false positive rate');
ylabel('true positive rate'); 
axis square; 

% separation
subplot(2, 2, 2);
l = result.averageRun.lambda; % p = result.averageRun.prototypes;
[lV, lD] = eig(l);
w = lV(:, end); % leading eigenvector of lambda
x = w' * X;
xWeights = x;
yMin = min(x);
yMax = max(x);
nSubjects = size(X, 2);
fill([1, 1, N1, N1], [yMin, yMax, yMax, yMin], [1.0, 0.9, 0.9], 'LineStyle', 'none'); hold on;
line([1, nSubjects], [0, 0], 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2); hold on;
plot(x, 'o-'); hold on; 
box off; xlim([1, size(X, 2)]); ylim([yMin, yMax]);
xlabel('patients / controls'); ylabel('w^T * X');
accuracy_n = mean((x(:) < 0) == y_TR(:));
accuracy_p = mean((x(:) > 0) == y_TR(:));
title(sprintf('accuracy = %1.2f', max(accuracy_n, accuracy_p)), 'FontWeight', 'normal');

% projections
subplot(2, 2, [3, 4]);
l = result.averageRun.lambda;
[lV, lD] = eig(l);
w = lV(:, end); % leading eigenvector of lambda
x = abs(V * S(1:K, 1:K) * w); % abs(V * w)
[~, ix] = sort(x);
x = x(ix);
plot(x, '.'); hold on; box off; xlim([0, length(V)]);
xlabel('metabolites'); ylabel('V * w');

% plot | GMLVQ default
% plot(result);

% save metabolites
bSaveMetabolites = 0;
if bSaveMetabolites == 1
  tCodes = names(ix);
  tCodes = tCodes(end:-1:1);
  xProjections = x;
  xProjections = xProjections(end:-1:1);
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_GMLVQ\\metabolites\\%s_names.txt', aFile);
  hFile = fopen(aFilename, 'wt');
  for i = 1:length(tCodes)
    fprintf(hFile, '%1.4f\t%s\n', xProjections(i), tCodes{i});
  end
  fclose(hFile);
end

% save weights
bSaveWeights = 0;
if bSaveWeights == 1
  labels = Y_labels;
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_GMLVQ\\weights\\%s_weights.mat', aFile);
  save(aFilename, 'xWeights', 'subjects', 'labels');
end

end % end

%-------------------------------------------------------------------------------