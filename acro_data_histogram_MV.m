%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_histogram_MV()

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
  x = isnan(data);

  xPatients = mean(x(:, pLabelsClass == 1), 2);
  xControls = mean(x(:, pLabelsClass == 0), 2);

  uPatients = unique(xPatients);
  uControls = unique(xControls);

  p = zeros(length(uPatients), length(uControls));
  for i = 1:length(uPatients)
    for j = 1:length(uControls)
      p(i, j) = sum(xPatients == uPatients(i) & xControls == uControls(j));
    end 
  end

  b = 0:(1.001 * mean(diff(uPatients))):1;
  hPatients = histc(xPatients, b);
  hControls = histc(xControls, b);

  % plot
  subplot(2, 2, 1 + (iRAW - 1) * 2); imagesc(uPatients, uControls, log10(p)); colorbar; 
  xlabel('Metabolites (controls)'); ylabel('Metabolites (patients)'); set(gca, 'YDir', 'normal')
  if bRAW == 1
    title('Fraction of missing values (original)', 'FontWeight', 'normal');
  else
    title('Fraction of missing values (with exclusion)', 'FontWeight', 'normal');
  end
  
  subplot(2, 2, 2 + (iRAW - 1) * 2); bar(b, [hPatients, hControls]); box off; legend({'Patients', 'Controls'});
  xlabel('Fraction of missing values'); ylabel('metabolites');
  if bRAW == 1
    title('Fraction of missing values (original)', 'FontWeight', 'normal');
  else
    title('Fraction of missing values (with exclusion)', 'FontWeight', 'normal');
  end
end

% save figure
aFilename = [aPath, '\\', 'histogram', '\\', aFile, '_histogram.png'];
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

end % end

%-------------------------------------------------------------------------------