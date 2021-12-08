%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function acro_data_histogram_feats()

% histogram of features

clc;

aPath = 'd:\\data\\acromegaly';
aFile = 'Lipids_P'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'
% note: biomarker 'HILIC_P' is detected ('M235T108')

% parameters
nMaxMissingValuesBySubjects = 0.30;
nMaxMissingValuesByFeatures = 0.50; 

% open figure
hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 1.75, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

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

  % nan mean
  x = data;
  xPatients = nanmean(x(:, pLabelsClass == 1), 2);
  xControls = nanmean(x(:, pLabelsClass == 0), 2);
  bLOG10 = 1;
  if bLOG10 == 1
    xPatients = log10(xPatients);
    xControls = log10(xControls);
  end
  p_val = ranksum(xPatients, xControls);

  xMin = min([xPatients(:); xControls(:)]);
  xMax = max([xPatients(:); xControls(:)]);

  nBins = 150;
  b = xMin:((xMax - xMin) / nBins):xMax;
  hPatients = histc(xPatients, b);
  hControls = histc(xControls, b);
  
  [~, i] = sort(xControls); 
  yPatients = xPatients(i);
  yControls = xControls(i);
  
  % plot
  subplot(2, 2, 1 + (iRAW - 1) * 2); plot(b, log10([hPatients, hControls]), 'x'); box off; legend({'Patients', 'Controls'});
  xlabel('log10(peak_intensities_of_metabolites)', 'Interpreter', 'none'); ylabel('log10(counts)'); 
  if bRAW == 1
    title(sprintf('Peak intensities (original) | p=%1.3f', p_val), 'FontWeight', 'normal');
  else
    title(sprintf('Peak intensities (with exclusion) | p=%1.3f', p_val), 'FontWeight', 'normal');
  end
  
  subplot(2, 2, 2 + (iRAW - 1) * 2); plot(yPatients, '-*'); hold on; plot(yControls, 'LineWidth', 2); box off;
  xlabel('Metabolites'); ylabel('log10(peak_intensities)', 'Interpreter', 'none'); xlim([1, length(yPatients)]);
  if bRAW == 1
    title('Sorted peak intensities (original)', 'FontWeight', 'normal');
  else
    title('Sorted peak intensities (with exclusion)', 'FontWeight', 'normal');
  end
end

% save figure
aFilename = [aPath, '\\', '_analysis', '\\', 'histogram', '\\', 'features', '\\', aFile, '_histogram.png'];
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

end % end

%-------------------------------------------------------------------------------