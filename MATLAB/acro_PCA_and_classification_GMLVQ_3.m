%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_PCA_and_classification_GMLVQ_3()

clc;

% add toolbox
addpath('gmlvq-v3-1');

aFile = 'HILIC_N'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

aClassifier = 'GMLVQ'; % 'SVM', 'GMLVQ'

bZScore = 1;

nTrainTestSplit = 0.75; % change to 75/25
% nConfidenceLimit = 0.95;

nRandomSamples = 100; % 5000

nPercentageOfSamples = 1.0; % nPercentageOfSamples *  min(patients), nPercentageOfSamples * controls
nEqualiseSamples = 1;

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

% svd
[~, S, V] = svd(cov(Y')); 
d = diag(S);

% PCA MATLAB | https://uk.mathworks.com/help/stats/pca.html
% [eigenvectors, scores, latent, tsquared, explained, mu] = pca(Y');

% reconstruct data
% PCA_reconstruct = scores * eigenvectors + mu;

% debug data reconstruction
% PCA = V' * Y;
% Y_R = V * PCA;
% debug^

% debug
% x = V' * Y; % j = 2; plot(x(j, :)); hold on; plot(scores(:, j));
% debug^

% explained variance
ev = cumsum(d) / sum(d);

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

% init
nMetabolites = size(V, 1);
pAccuracies = zeros(nRandomSamples, 2);
pWeights = zeros(K, nRandomSamples, 2);
pProducts = zeros(nMetabolites, nRandomSamples, 2);
% loop shuffles
for bShuffle = 0:1
  % init splits
  nSplit = nTrainTestSplit;
  fprintf(1, 'split: %d/%d | shuffled (%d)\n', round(100 * nSplit), round(100 * (1 - nSplit)), bShuffle);
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
  
  % random sampling
  for iRandomSample = 1:nRandomSamples
    fprintf(1, '%d / %d\n', iRandomSample, nRandomSamples);
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
    if bShuffle == 1
      y_TR = y_TR(randperm(length(y_TR)));
      y_TS = y_TS(randperm(length(y_TS)));
    end

    % z-score
    if bZScore == 1
      X_TR = (X_TR - mean(X_TR, 2)) ./ std(X_TR, 0, 2);
      X_TS = (X_TS - mean(X_TS, 2)) ./ std(X_TS, 0, 2);
    end
    
    % classification GMLVQ
    if strcmp(aClassifier, 'GMLVQ')
      % single run
      bSingleRun = 0;
      if bSingleRun == 1
        GMLVQ_model = GMLVQ.GMLVQ(X_TR', (y_TR + 1)', GMLVQ.Parameters(), 30);
        result = GMLVQ_model.runSingle();
        validationData = GMLVQ.DataPair(X_TS', (y_TS + 1)');
        [~, GMLVQ_labels, ~, ~] = result.run.classify(validationData);
        GMLVQ_accuracy = mean(GMLVQ_labels(:) == (y_TS(:) + 1));
        lambda = result.run.lambda;
        prototypes = result.run.prototypes;
      else
        % validation
        GMLVQ_model = GMLVQ.GMLVQ([X_TR, X_TS]', ([y_TR, y_TS] + 1)', GMLVQ.Parameters(), 30); % 30 == totalsteps
        result = GMLVQ_model.runValidation(10, 25); % 10 validation runs, 25% of examples per class
        lambda = result.averageRun.lambda;
        prototypes = result.averageRun.prototypes;
      end

      [e_vec, e_val] = eig(lambda);  

      % plot(result);

      % init
      proj = e_vec(:, end); % leading eigenvector of lambda
      proj = proj / vecnorm(proj);
      pAccuracies(iRandomSample, bShuffle + 1) = GMLVQ_accuracy;
      pWeights(:, iRandomSample, bShuffle + 1) = proj;
      pProducts(:, iRandomSample, bShuffle + 1) = V * S(1:K, 1:K) * proj; 

    elseif strcmp(aClassifier, 'SVM')
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
      pAccuracies(iRandomSample, bShuffle + 1) = SVM_accuracy;
      pWeights(:, iRandomSample, bShuffle + 1) = beta; 
      pProducts(:, iRandomSample, bShuffle + 1) = V * S(1:K, 1:K) * beta; 
    else
      fprintf(1, 'ERROR: classifier ''%s'' is not implemented.\n', aClassifier);
      return
    end
  end
  % histogram
  pBins = 0.0:0.1:1;
  subplot(2, 3, 1);
  h = histc(pAccuracies(:, bShuffle + 1), pBins); hold on;
  plot(pBins, h / nRandomSamples, 'o-'); xlim([0.25, 0.95]); box off; 
  xlabel('accuracy'); ylabel('counts / sum');
  title(sprintf('hist | train/test: %d/%d', 100 * nSplit, 100 * (1 - nSplit)), 'FontWeight', 'normal');

  % linear SVM
  if strcmp(aClassifier, 'SVM')
    % weights
    subplot(2, 3, 2);
    w = squeeze(mean(pWeights(:, :, bShuffle + 1), 2));
    plot(w, 'o-'); hold on; 
    
    xlim([0.5, K + 0.5]); box off; 
    xlabel('PC index'); ylabel('weights');
    title(sprintf('w | K=%d | train/test: %d/%d', K, 100 * nSplit, 100 * (1 - nSplit)), 'FontWeight', 'normal');

    % separation
    w = squeeze(mean(pWeights(:, :, bShuffle + 1), 2));
    x = w' * X;
    yMin = min(x);
    yMax = max(x);
    subplot(2, 3, 3);
    if bShuffle == 0
      fill([1, 1, N1, N1], [yMin, yMax, yMax, yMin], [1.0, 0.9, 0.9], 'LineStyle', 'none'); hold on;
    end
    plot(x, 'o-'); hold on; 
    if bShuffle == 0
      box off; xlim([1, size(X, 2)]); ylim([yMin, yMax]);
      xlabel('patients / controls'); ylabel('w^T * X');
    end
    title(sprintf('w^T * X | K=%d | train/test: %d/%d', K, 100 * nSplit, 100 * (1 - nSplit)), 'FontWeight', 'normal');
    
    % products
    subplot(2, 3, [4, 6]);

    w = squeeze(mean(pWeights(:, :, bShuffle + 1), 2));
    x = abs(V * S(1:K, 1:K) * w); % abs(V * w)
    if bShuffle == 0
      [~, ix] = sort(x);
    end
    x = x(ix);
  
    gx = x(end:-1:1);
    cx = cumsum(gx) / sum(gx);
    gx = gx(cx < 0.10);
  
    % min distance
    bMinDistance = 0;
    if bMinDistance == 1
      yx = gx / max(gx);
      xx = (1:length(x))' / length(x); 
      dx = sqrt(xx .^ 2 + yx .^ 2);
      plot(yx); hold on; plot(dx);
    end
  
    % fit 
    % fx = x(end:-1:1);
    % fL = 500;
    % fN = length(x) - fL;
    % fd = zeros(fN, 1);
    % for i = 1:fN
    %   iBeg = i;
    %   iEnd = iBeg + fL - 1;
    %   gx = fx(iBeg:iEnd);
    %   [u, gof] = support_fit_spectra((1:fL)', gx, 'lin');
    %   fd(i) = gof;
    % end
  
    % powerlaw fit
    % [u, gof] = support_fit_spectra((1:length(gx))', gx, 'pow');
    % plot(diff(u) / sum(diff(u)));
  
    % title
    plot(x); hold on; box off; xlim([0, length(V)]);
    title(sprintf('V * w | train/test: %d/%d', 100 * nSplit, 100 * (1 - nSplit)), 'FontWeight', 'normal');
    xlabel('metabolites'); ylabel('V * w');
  end

  % GMLVQ
  if strcmp(aClassifier, 'GMLVQ')
    % weights
    subplot(2, 3, 2);
    w = squeeze(mean(pWeights(:, :, bShuffle + 1), 2));
    plot(w, 'o-'); hold on; 
    nSignificants = [];
    xlim([0.5, K + 0.5]); box off; 
    xlabel('PC index'); ylabel('weights');
    if ~isempty(nSignificants)
      title(sprintf('w | K=%d | train/test: %d/%d | %d', K, 100 * nSplit, 100 * (1 - nSplit), nSignificants), 'FontWeight', 'normal');
    else
      title(sprintf('w | K=%d | train/test: %d/%d', K, 100 * nSplit, 100 * (1 - nSplit)), 'FontWeight', 'normal');
    end
    % separation
    w = squeeze(mean(pWeights(:, :, bShuffle + 1), 2));
    x = w' * X;
    yMin = min(x);
    yMax = max(x);
    subplot(2, 3, 3);
    if bShuffle == 0
      fill([1, 1, N1, N1], [yMin, yMax, yMax, yMin], [1.0, 0.9, 0.9], 'LineStyle', 'none'); hold on;
    end
    plot(x, 'o-'); hold on; 
    if bShuffle == 0
      box off; xlim([1, size(X, 2)]); ylim([yMin, yMax]);
      xlabel('patients / controls'); ylabel('w^T * X');
    end
    title(sprintf('w^T * X | K=%d | train/test: %d/%d', K, 100 * nSplit, 100 * (1 - nSplit)), 'FontWeight', 'normal');

    % products CI
    subplot(2, 3, [4, 6]);
    % products
    w = squeeze(mean(pWeights(:, :, bShuffle + 1), 2));
    x = abs(V * S(1:K, 1:K) * w); % abs(V * w)
    if bShuffle == 0
      [~, ix] = sort(x);
    end
    x = x(ix);
    plot(x); hold on; box off; xlim([0, length(V)]);
    nSignificants = [];
    % title
    if ~isempty(nSignificants)
      title(sprintf('V * w | train/test: %d/%d | %d', 100 * nSplit, 100 * (1 - nSplit), nSignificants), 'FontWeight', 'normal');
    else
      title(sprintf('V * w | train/test: %d/%d', 100 * nSplit, 100 * (1 - nSplit)), 'FontWeight', 'normal');
    end
    xlabel('metabolites'); ylabel('V * w');
  end
end

fprintf(1, 'K = %d\n', K);

% save weights
bSaveWeights = 0;
if bSaveWeights == 1
  w = squeeze(mean(pWeights(:, :, 1), 2));
  x = w' * X;
  labels = Y_labels;
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\weights\\%s_%d%d.mat', aFile, 100 * nSplit, 100 * (1 - nSplit));
  save(aFilename, 'w', 'x', 'subjects', 'labels');
end

end % end

%-------------------------------------------------------------------------------