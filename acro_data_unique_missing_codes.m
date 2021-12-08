%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_unique_missing_codes()

% histogram of missing values

clc;

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
  bRAW = 2 - iRAW;
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

  % missing values to +/-1
  x = double(isnan(data));
  
  nFeatures = size(x, 1);
  nSubjects = size(x, 2);
  pCodes = zeros(nFeatures, 1);
  tCodes = cell(nFeatures, 1);
  W = 2 .^ (0:(nSubjects - 1));
  for iFeature = 1:nFeatures
    pCodes(iFeature) = sum(x(iFeature, :) .* W);
    tCodes{iFeature} = strrep(num2str(x(iFeature, :)), ' ', '');
  end

  pUniqueCodes = unique(pCodes);
  tUniqueCodes = unique(tCodes);
  if pUniqueCodes(1) == 0
    pUniqueCodes(1) = [];
    tUniqueCodes(1) = [];
  end
  nUniqueCodes = length(pUniqueCodes);
  pCodesCounts = zeros(nUniqueCodes, 1);
  for i = 1:nUniqueCodes
    pCodesCounts(i) = sum(pCodes == pUniqueCodes(i));
  end
  
  nMinRecurrence = 5;
  pUniqueCodesRepeated = pUniqueCodes(pCodesCounts > nMinRecurrence);
  tUniqueCodesRepeated = tUniqueCodes(pCodesCounts > nMinRecurrence);
  pUniqueCodesRepeatedCounts = pCodesCounts(pCodesCounts > nMinRecurrence);
  
  [~, i] = sort(pUniqueCodesRepeatedCounts);
  pUniqueCodesRepeated = pUniqueCodesRepeated(i);
  tUniqueCodesRepeated = tUniqueCodesRepeated(i);
  pUniqueCodesRepeatedCounts = pUniqueCodesRepeatedCounts(i);
  
  nUniqueFeatures = length(pUniqueCodesRepeatedCounts);
  X = zeros(nUniqueFeatures, nSubjects);
  for iFeature = 1:nUniqueFeatures
    for i = 1:nSubjects
      X(iFeature, i) = str2double(tUniqueCodesRepeated{iFeature}(i));
    end
  end
  
  % log10
  pUniqueCodesRepeatedCounts = log10(pUniqueCodesRepeatedCounts);
  
  % plot
  subplot(2, 4, [1, 3] + (iRAW - 1) * 4); imagesc(X, [-2, 2]); colormap('jet'); box off; colorbar;
  xlabel('Subjects'); ylabel('Unique patterns');
  subplot(2, 4, 4 + (iRAW - 1) * 4); imagesc(pUniqueCodesRepeatedCounts, [0, max(pUniqueCodesRepeatedCounts)]); hcb = colorbar;
  hcb.Label.String = 'log10(recurrence)';
  ylabel('Unique patterns');
end

% save figure
aFilename = [aPath, '\\', '_analysis', '\\', 'unique_codes', '\\', aFile, '_codes.png'];
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

end % end

%-------------------------------------------------------------------------------