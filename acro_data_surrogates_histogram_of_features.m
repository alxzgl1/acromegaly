%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_surrogates_histogram_of_features()

clc;

aPath = 'd:\\data\\acromegaly';
aFile = 'Lipids_N'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

aSurrogates = 'by_features'; % 'by_subjects', 'by_features'

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

% missing values 
x = ~isnan(data);

% randomise subject order
nSubjects = size(x, 2);
nFeatures = size(x, 1);
nSurrogates = 100;
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

b = unique([sum(x, 2) / nSubjects; sum(y, 2) / nSubjects]);
hx = histc(sum(x, 2) / nSubjects, b);
hy = histc(sum(y, 2) / nSubjects, b);

HY = zeros(length(b), nSurrogates);
for iSurrogate = 1:nSurrogates
  v = Y(:, :, iSurrogate);
  hv = histc(sum(v, 2) / nSubjects, b);
  HY(:, iSurrogate) = hv; 
end

% open figure
hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 1.5, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 
% plot
subplot(2, 3, 1); imagesc(x, [0, 1]); set(gca, 'YDir', 'normal'); colorbar;
xlabel('subjects'); ylabel('metabolites'); title('A1 | Original data (missing is 0)', 'FontWeight', 'normal');
subplot(2, 3, 4); imagesc(y, [0, 1]); set(gca, 'YDir', 'normal'); colorbar;
xlabel('subjects'); ylabel('metabolites'); title(['A2 | Shuffled ', aSurrogates, ' (1 run)'], 'FontWeight', 'normal', 'Interpreter', 'none');

% plot
subplot(2, 3, 2); plot(sum(x, 2) / nSubjects, 1:nFeatures); box off; ylim([1, nFeatures]);
xlabel('Percentage of non-missing values'); ylabel('metabolites'); title('B1 | Percentage of non-missing values', 'FontWeight', 'normal');
subplot(2, 3, 5); plot(sum(y, 2) / nSubjects, 1:nFeatures); box off; ylim([1, nFeatures]);
xlabel('Percentage of non-missing values'); ylabel('metabolites'); title(['B2 | Shuffled ', aSurrogates, ' (1 run)'], 'FontWeight', 'normal', 'Interpreter', 'none');

% histogram
yMin = min([hx(:); HY(:)]);
yMax = max([hx(:); HY(:)]);
subplot(2, 3, 3); bar(b, hx); ylim([yMin, yMax]); box off;
xlabel('Percentage of non-missing values'); ylabel('Counts'); title('C1 | Histogram of non-missing values B1', 'FontWeight', 'normal');
for iSurrogate = 1:nSurrogates
  subplot(2, 3, 6); bar(b, HY(:, iSurrogate)); hold on; 
end
ylim([yMin, yMax]); box off;
xlabel('Percentage of non-missing values'); ylabel('Counts'); title('C2 | Histogram of non-missing values B2 (100 runs)', 'FontWeight', 'normal');
% save figure
aFile = sprintf('%s_MV%d%d', aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
aFilename = [aPath, '\\', '_analysis', '\\', 'histogram_of_features', '\\', 'surrogates_', aSurrogates, '\\', aFile, '.png'];
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

end % end

%-------------------------------------------------------------------------------