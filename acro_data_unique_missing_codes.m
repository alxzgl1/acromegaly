%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_unique_missing_codes()

clc;

nMinCodeRecurrence = 5; % 5 (default)

bShuffle = 0; % 0 (default), if 1, run again until NO error, and try reducing nMinCodeRecurrence

aPath = 'd:\\data\\acromegaly';
aFile = 'Lipids_P'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'
% note: biomarker 'HILIC_P' is detected ('M235T108')

% parameters
nMaxMissingValuesBySubjects = 0.30;
nMaxMissingValuesByFeatures = 0.50; 

% open figure
hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

% loop
for iRAW = 1:2
  bRAW = 2 - iRAW; if bShuffle == 1, bRAW = 0; end
  % load data
  if bRAW == 1
    aFilename = sprintf('%s\\import\\%s_data.mat', aPath, aFile);
    load(aFilename, 'data');
    aFilename = sprintf('%s\\import\\%s_names.mat', aPath, aFile);
    load(aFilename, 'names');
    aFilename = sprintf('%s\\import\\%s_labels.mat', aPath, aFile);
    load(aFilename, 'labels');
  else
    aFilename = sprintf('%s\\import\\%s_data_MV%d%d.mat', aPath, aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
    load(aFilename, 'data');
    aFilename = sprintf('%s\\import\\%s_names_MV%d%d.mat', aPath, aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
    load(aFilename, 'names');
    aFilename = sprintf('%s\\import\\%s_labels_MV%d%d.mat', aPath, aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
    load(aFilename, 'labels');
  end

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
  x = double(~isnan(data));

  % shuffle
  nFeatures = size(x, 1);
  nSubjects = size(x, 2);
  if bShuffle == 1
    if iRAW == 2
      for i = 1:nSubjects
        x(:, i) = x(randperm(nFeatures), i);
      end
    end
  end
  
  nFeatures = size(x, 1);
  nSubjects = size(x, 2);
  tCodes = cell(nFeatures, 1);
  for iFeature = 1:nFeatures
    tCodes{iFeature} = strrep(num2str(x(iFeature, :)), ' ', '');
  end

  % unique codes
  tUniqueCodes = unique(tCodes);
  bFeaturesFull = 0;
  if sum(double(tUniqueCodes{end}) / double('1')) == nSubjects % non-missing
    bFeaturesFull = 1; % tUniqueCodes(end) = [];
  end
  nUniqueCodes = length(tUniqueCodes);
  pCodesCounts = zeros(nUniqueCodes, 1);
  for i = 1:nUniqueCodes
    pCodesCounts(i) = sum(strcmp(tCodes, tUniqueCodes{i}));
  end
  
  tUniqueCodesRepeated = tUniqueCodes(pCodesCounts > nMinCodeRecurrence);
  pUniqueCodesRepeatedCounts = pCodesCounts(pCodesCounts > nMinCodeRecurrence);
  
  [~, i] = sort(pUniqueCodesRepeatedCounts);
  tUniqueCodesRepeated = tUniqueCodesRepeated(i);
  pUniqueCodesRepeatedCounts = pUniqueCodesRepeatedCounts(i);

  % exclude non-missing values
  if bFeaturesFull == 1
    nFeaturesFull = round(100 * (pUniqueCodesRepeatedCounts(end) / nFeatures));
    tUniqueCodesRepeated(end) = [];
    pUniqueCodesRepeatedCounts(end) = [];
    nFeaturesAboveMin = round(100 * (sum(pUniqueCodesRepeatedCounts) / nFeatures));
    nFeaturesBelowMin = 100 - nFeaturesFull - nFeaturesAboveMin;
  else
    nFeaturesFull = 0;
    nFeaturesAboveMin = round(100 * (sum(pUniqueCodesRepeatedCounts) / nFeatures));
    nFeaturesBelowMin = 100 - nFeaturesFull - nFeaturesAboveMin;
  end
  
  nUniqueFeatures = length(pUniqueCodesRepeatedCounts);
  X = zeros(nUniqueFeatures, nSubjects);
  for iFeature = 1:nUniqueFeatures
    for i = 1:nSubjects
      X(iFeature, i) = str2double(tUniqueCodesRepeated{iFeature}(i));
    end
  end
  
  % log10
  bLOG10 = 0;
  if bLOG10 == 1
    pUniqueCodesRepeatedCounts = log10(pUniqueCodesRepeatedCounts);
  end
  
  % plot
  subplot(2, 4, [1, 3] + (iRAW - 1) * 4); imagesc(X, [0, 1]); box off; colorbar; 
  xlabel('Subjects'); ylabel('Unique patterns (missing=blue)'); set(gca, 'YDir', 'normal');
  if bShuffle == 1
    if iRAW == 1
      title(sprintf('Codes after exclusion | nonmissing=%d%%, above min=%d%%, below min=%d%%', nFeaturesFull, nFeaturesAboveMin, nFeaturesBelowMin), 'FontWeight', 'normal');
    else
      title(sprintf('Codes after exclusion shuffled | nonmissing=%d%%, above min=%d%%, below min=%d%%', nFeaturesFull, nFeaturesAboveMin, nFeaturesBelowMin), 'FontWeight', 'normal');
    end
  else
    if bRAW == 1
      title(sprintf('Codes before exclusion | nonmissing=%d%%, above min=%d%%, below min=%d%%', nFeaturesFull, nFeaturesAboveMin, nFeaturesBelowMin), 'FontWeight', 'normal');
    else
      title(sprintf('Codes after exclusion | nonmissing=%d%%, above min=%d%%, below min=%d%%', nFeaturesFull, nFeaturesAboveMin, nFeaturesBelowMin), 'FontWeight', 'normal');
    end
  end
  % add grid
  rows = size(X, 1) + 2;
  columns = size(X, 2) + 2;
  for row = 1:rows
    line([0, columns + 1], [row - 0.5, row - 0.5], 'Color', 'k', 'LineWidth', 0.5);
  end
  for col = 1:columns
    line([col - 0.5, col - 0.5], [0, rows + 1], 'Color', 'k', 'LineWidth', 0.5);
  end
  % bar
  subplot(2, 4, 4 + (iRAW - 1) * 4); barh(1:length(pUniqueCodesRepeatedCounts), pUniqueCodesRepeatedCounts); 
  ylim([0.5, length(pUniqueCodesRepeatedCounts) + 0.5]); box off;
  ylabel('Unique patterns'); 
  title(sprintf('Min recurrence is %d', nMinCodeRecurrence), 'FontWeight', 'normal');
  if bLOG10 == 1
    xlabel('log10(recurrence)'); 
  else
    xlabel('Recurrence counts'); 
  end
end

% save figure
if bShuffle == 1
  aFilename = [aPath, '\\', '_analysis', '\\', 'unique_codes', '\\', 'shuffled', '\\', aFile, '_codes_shuffled.png'];
else
  aFilename = [aPath, '\\', '_analysis', '\\', 'unique_codes', '\\', aFile, '_codes.png'];
end
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

end % end

%-------------------------------------------------------------------------------