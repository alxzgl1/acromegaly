%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_surrogates()

clc;

aPath = 'd:\\data\\acromegaly';
aFile = 'Lipids_P'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

% note: biomarker 'HILIC_P' is detected ('M235T108')

% parameters
nMaxMissingValuesBySubjects = 0.30;
nMaxMissingValuesByFeatures = 0.50; 

% load data and header
aFilename = sprintf('%s\\import\\%s_data_MV%d%d.mat', aPath, aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
load(aFilename, 'data');
aFilename = sprintf('%s\\import\\%s_names_MV%d%d.mat', aPath, aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
load(aFilename, 'names');
aFilename = sprintf('%s\\import\\%s_labels_MV%d%d.mat', aPath, aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
load(aFilename, 'labels');

% CHECK
if size(data, 1) ~= size(names, 1)
  fprintf(1, 'WARNING: mismatch between data (columns) and names.\n');
  return
end
if size(data, 2) ~= size(labels, 2) - 1
  fprintf(1, 'WARNING: mismatch between data (rows) and labels.\n');
  return
end

% parse labels
tLabelsID = labels(1, 2:end);
tLabelsClass = labels(contains(labels(:, 1), 'Class'), 2:end);
pLabelsClass = contains(tLabelsClass, 'Acromegaly');
tLabelsSex = labels(contains(labels(:, 1), 'Sex'), 2:end);
tLabelsAge = labels(contains(labels(:, 1), 'Age'), 2:end);
tLabelsBMI = labels(contains(labels(:, 1), 'BMI'), 2:end);

% CHECK
if length(unique(tLabelsClass)) ~= 2
  fprintf(1, 'WARNING: number of classes should exactly 2 (not %d as detected).\n', length(unique(tLabelsClass)));
  return
end

% missing values to +/-1
x = sign(0.5 - isnan(data));

% remove features with non-missing values
p = mean(x, 2);
x(p == 1, :) = [];
% names(p == 1) = [];

% randomise subject order
nSubjects = size(x, 2);
nFeatures = size(x, 1);
nSurrogates = 100;
Y = zeros(nFeatures, nSubjects, nSurrogates);
for iSurrogate = 1:nSurrogates
  for i = 1:nFeatures
    Y(i, :, iSurrogate) = x(i, randperm(nSubjects));
  end
end
y = Y(:, :, 1); % 1 run

% randomise feature order
s = x(randperm(nFeatures), :);

cx = corr(x'); 
cs = corr(s'); 
cy = corr(y'); 

PY = zeros(nFeatures, nSurrogates);
for iSurrogate = 1:nSurrogates
  v = Y(:, :, iSurrogate);
  PY(:, iSurrogate) = sort(sum(corr(v'), 2)); 
end

% layout 3x3
bLayout3x3 = 0;
if bLayout3x3 == 1 
  % open figure
  hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 1.5, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 
  % plot
  subplot(3, 3, 1); imagesc(x, [-1, 1]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('subjects'); ylabel('metabolites'); title('A1 | Original data (missing (-1), existing (+1))', 'FontWeight', 'normal');
  subplot(3, 3, 2); imagesc(s, [-1, 1]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('subjects'); ylabel('metabolites'); title('A2 | Shuffled order of metabolites', 'FontWeight', 'normal');
  subplot(3, 3, 3); imagesc(y, [-1, 1]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('subjects'); ylabel('metabolites'); title('A3 | Shuffled order of subjects (1 run)', 'FontWeight', 'normal');
  % correlation
  subplot(3, 3, 4); imagesc(cx, [-0.3, 0.3]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('metabolites'); ylabel('metabolites'); title('B1 | Correlation matrix of A1', 'FontWeight', 'normal');
  subplot(3, 3, 5); imagesc(cs, [-0.3, 0.3]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('metabolites'); ylabel('metabolites'); title('B2 | Correlation matrix of A2', 'FontWeight', 'normal');
  subplot(3, 3, 6); imagesc(cy, [-0.3, 0.3]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('metabolites'); ylabel('metabolites'); title('B3 | Correlation matrix of A3 (1 run)', 'FontWeight', 'normal');
  % sorted sum
  px = sort(sum(cx, 2));
  ps = sort(sum(cs, 2));
  xMin = min([px(:); ps(:); PY(:)]);
  xMax = max([px(:); ps(:); PY(:)]);
  subplot(3, 3, 7); plot(px, 1:nFeatures, 'k.'); xlim([xMin, xMax]);
  xlabel('sorted sum'); ylabel('metabolites'); title('C1 | Sorted sum of Correlation matrix B1', 'FontWeight', 'normal');
  subplot(3, 3, 8); plot(ps, 1:nFeatures, 'k.'); xlim([xMin, xMax]);
  xlabel('sorted sum'); ylabel('metabolites'); title('C2 | Sorted sum of Correlation matrix B2', 'FontWeight', 'normal');
  for iSurrogate = 1:nSurrogates
    subplot(3, 3, 9); plot(PY(:, iSurrogate), 1:nFeatures, 'k.'); hold on; xlim([xMin, xMax]);
  end
  xlabel('sorted sum'); ylabel('metabolites'); title('C3 | Sorted sum of Correlation matrix B3 (100 runs)', 'FontWeight', 'normal');
  % save figure
  aFile = sprintf('%s_MV%d%d', aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
  aFilename = [aPath, '\\', '_analysis', '\\', 'surrogates', '\\', aFile, '.png'];
  print(hFigure, aFilename, '-dpng', '-r300');
  close(hFigure);
end

% layout 2x3
bLayout2x3 = 1;
if bLayout2x3 == 1 
  % open figure
  hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 1.5, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 
  % plot
  subplot(2, 3, 1); imagesc(x, [-1, 1]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('subjects'); ylabel('metabolites'); title('A1 | Original data (missing/existing)', 'FontWeight', 'normal');
  subplot(2, 3, 4); imagesc(y, [-1, 1]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('subjects'); ylabel('metabolites'); title('A3 | Shuffled order of subjects (1 run)', 'FontWeight', 'normal');
  % correlation
  subplot(2, 3, 2); imagesc(cx, [-0.2, 0.2]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('metabolites'); ylabel('metabolites'); title('B1 | Correlation matrix of A1', 'FontWeight', 'normal');
  subplot(2, 3, 5); imagesc(cy, [-0.2, 0.2]); set(gca, 'YDir', 'normal'); colorbar;
  xlabel('metabolites'); ylabel('metabolites'); title('B3 | Correlation matrix of A3 (1 run)', 'FontWeight', 'normal');
  % sorted sum
  px = sort(sum(cx, 2));
  ps = sort(sum(cs, 2));
  xMin = min([px(:); ps(:); PY(:)]);
  xMax = max([px(:); ps(:); PY(:)]);
  subplot(2, 3, 3); plot(px, 1:nFeatures, 'k.'); xlim([xMin, xMax]);
  xlabel('sorted sum'); ylabel('metabolites'); title('C1 | Sorted sum of Correlation matrix B1', 'FontWeight', 'normal');
  for iSurrogate = 1:nSurrogates
    subplot(2, 3, 6); plot(PY(:, iSurrogate), 1:nFeatures, 'k.'); hold on; xlim([xMin, xMax]);
  end
  xlabel('sorted sum'); ylabel('metabolites'); title('C3 | Sorted sum of Correlation matrix B3 (100 runs)', 'FontWeight', 'normal');
  % save figure
  aFile = sprintf('%s_MV%d%d', aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
  aFilename = [aPath, '\\', '_analysis', '\\', 'surrogates', '\\', aFile, '.png'];
  print(hFigure, aFilename, '-dpng', '-r300');
  close(hFigure);
end

end % end

%-------------------------------------------------------------------------------