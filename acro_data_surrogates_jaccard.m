%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_surrogates_jaccard()

clc;

aPath = 'd:\\data\\acromegaly';
aFile = 'Lipids_P'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

aSurrogates = 'by_features'; % 'by_subjects', 'by_features'

bDivideBySum = 1; % 1 (default)

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
% p = mean(x, 2);
% x(p == 1, :) = [];
% names(p == 1) = [];

% randomise subject order
nSubjects = size(x, 2);
nFeatures = size(x, 1);
nSurrogates = 5;
Y = zeros(nFeatures, nSubjects, nSurrogates);
if strcmp(aSurrogates, 'by_subjects')
  for iSurrogate = 1:nSurrogates
    for i = 1:nFeatures
      Y(i, :, iSurrogate) = x(i, randperm(nSubjects));
    end
  end
end
if strcmp(aSurrogates, 'by_features')
  for iSurrogate = 1:nSurrogates
    for i = 1:nSubjects
      Y(:, i, iSurrogate) = x(randperm(nFeatures), i);
    end
  end
end
y = Y(:, :, 1); % 1 run

cx = jaccard(x', bDivideBySum); cx = cx .* (1 - eye(nFeatures)); 
cy = jaccard(y', bDivideBySum); cy = cy .* (1 - eye(nFeatures)); 

nBins = 100;
pcx = cx(:); 
pcy = cy(:); 
xMin = min([pcx(:); pcy(:)]);
xMax = max([pcx(:); pcy(:)]);
b = xMin:((xMax - xMin) / nBins):xMax;
hx = histc(pcx, b);
hy = histc(pcy, b);

HY = zeros(length(b), nSurrogates);
for iSurrogate = 1:nSurrogates
  v = Y(:, :, iSurrogate);
  cv = jaccard(v', bDivideBySum); cv = cv .* (1 - eye(nFeatures)); 
  pcv = cv(:); 
  hv = histc(pcv, b);
  HY(:, iSurrogate) = hv; 
end

% open figure
hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 1.5, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 
% plot
subplot(2, 3, 1); imagesc(x, [0, 1]); set(gca, 'YDir', 'normal'); colorbar;
xlabel('subjects'); ylabel('metabolites'); title('A1 | Original data (missing/existing)', 'FontWeight', 'normal');
subplot(2, 3, 4); imagesc(y, [0, 1]); set(gca, 'YDir', 'normal'); colorbar;
xlabel('subjects'); ylabel('metabolites'); title(['A2 | Shuffled ', aSurrogates, ' (1 run)'], 'FontWeight', 'normal', 'Interpreter', 'none');
% correlation
subplot(2, 3, 2); imagesc(cx, [0, 1]); set(gca, 'YDir', 'normal'); colorbar;
xlabel('metabolites'); ylabel('metabolites'); title('B1 | Correlation matrix of A1', 'FontWeight', 'normal');
subplot(2, 3, 5); imagesc(cy, [0, 1]); set(gca, 'YDir', 'normal'); colorbar;
xlabel('metabolites'); ylabel('metabolites'); title('B2 | Correlation matrix of A2 (1 run)', 'FontWeight', 'normal');
% histogram
yMin = min([hx(:); HY(:)]);
yMax = max([hx(:); HY(:)]);
subplot(2, 3, 3); bar(b, hx); ylim([yMin, yMax]); box off;
xlabel('Correlation values'); ylabel('Counts'); title('C1 | Histogram of Correlation matrix B1', 'FontWeight', 'normal');
for iSurrogate = 1:nSurrogates
  subplot(2, 3, 6); bar(b, HY(:, iSurrogate)); hold on; 
end
ylim([yMin, yMax]); box off;
xlabel('Correlation values'); ylabel('Counts'); title('C2 | Histogram of Correlation matrix B2 (100 runs)', 'FontWeight', 'normal');
% save figure
aFile = sprintf('%s_MV%d%d', aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
if bDivideBySum == 1, aFile = [aFile, '_SUM']; end
aFilename = [aPath, '\\', '_analysis', '\\', 'jaccard', '\\', 'surrogates_', aSurrogates, '\\', aFile, '.png'];
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

end % end

%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function J = jaccard(x, bDivideBySum)

x = x > 0;

nSubjects = size(x, 1);
nFeatures = size(x, 2);

J = zeros(nFeatures, nFeatures);
for i = 1:nFeatures
  for j = (i + 1):nFeatures
    if bDivideBySum == 1
      J(i, j) = sum(x(:, i) & x(:, j)) / nSubjects;
    else
      J(i, j) = sum(x(:, i) & x(:, j)) / sum(x(:, i) | x(:, j));
    end
  end
end
J = J + J';

end % end

%-------------------------------------------------------------------------------