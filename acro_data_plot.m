%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_plot()

clc;

aPath = 'd:\\data\\acromegaly';
aFile = 'HILIC_P'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

% parameters
nMaxMissingValuesBySubjects = 0.30;
nMaxMissingValuesByFeatures = 0.50; 
nMaxMissingValuesBiomarkers = 0.825; % keep above 80% separately as potential biomarker

% load data and header
aFilename = [aPath, '\\', 'import', '\\', aFile, '_data.mat'];
load(aFilename, 'data');
aFilename = [aPath, '\\', 'import', '\\', aFile, '_names.mat'];
load(aFilename, 'names');
aFilename = [aPath, '\\', 'import', '\\', aFile, '_labels.mat'];
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
tLabelsClass = labels(contains(labels(:, 1), 'Class'), 2:end);
pLabelsClass = contains(tLabelsClass, 'Acromegaly');

% CHECK
if length(unique(tLabelsClass)) ~= 2
  fprintf(1, 'WARNING: number of classes should exactly 2 (not %d as detected).\n', length(unique(tLabelsClass)));
  return
end

% keep data
data_x = data;
labels_x = pLabelsClass;
nMetabolites_x = length(labels_x);
         
% count missing values by subjects       
nFeatures = size(data, 1);
pMissingValuesBySubjects = sum(isnan(data)) / nFeatures;
pBadSubjects = pMissingValuesBySubjects > nMaxMissingValuesBySubjects;

nMaxPatients = sum(pLabelsClass == 1);
nMaxControls = sum(pLabelsClass == 0);

% exclude bad subjects 
data = data(:, pBadSubjects == 0);
pLabelsClass = pLabelsClass(pBadSubjects == 0);

nPatients = sum(pLabelsClass == 1);
nControls = sum(pLabelsClass == 0);

% status
fprintf(1, 'STATUS: %d%% (patients) and %d%% (controls) excluded.\n', round(100 * (1 - nPatients / nMaxPatients)), round(100 * (1 - nControls / nMaxControls)));
 
% count missing values by features
pMissingValuesByFeatures = [sum(isnan(data(:, pLabelsClass == 1)), 2) / nPatients, ...
                            sum(isnan(data(:, pLabelsClass == 0)), 2) / nControls];

pBadFeatures = pMissingValuesByFeatures(:, 1) > nMaxMissingValuesByFeatures | ...
               pMissingValuesByFeatures(:, 2) > nMaxMissingValuesByFeatures;

% check potential biomarkers
bCheckBiomarkers = 1;
if bCheckBiomarkers == 1
  pBiomarkers = ((pMissingValuesByFeatures(:, 1) < (1 - nMaxMissingValuesBiomarkers)) & (pMissingValuesByFeatures(:, 2) > nMaxMissingValuesBiomarkers)) | ...
                ((pMissingValuesByFeatures(:, 1) > nMaxMissingValuesBiomarkers) & (pMissingValuesByFeatures(:, 2) < (1 - nMaxMissingValuesBiomarkers)));
  if sum(pBiomarkers) > 0
    fprintf(1, 'WARNING: possible biomarkers detected.\n');
    p = find(pBiomarkers > 0);
    for i = 1:length(p)
      aFeature = names{p(i)};
      fprintf(1, '%s | missing values: %d%% (patients) vs %d%% (controls).\n', aFeature, ...
        round(100 * pMissingValuesByFeatures(p(i), 1)), round(100 * pMissingValuesByFeatures(p(i), 2)));
    end
  end
  % handle data
  xPatients = round(100 * pMissingValuesByFeatures(:, 1));
  xControls = round(100 * pMissingValuesByFeatures(:, 2));
  [~, i] = sort(xControls, 'descend');
  nMetabolites = length(i);
  xPatients = xPatients(i);
  xControls = xControls(i);
  tMetabolites = names(i);
  pBiomarkers = ((xPatients < 100 * (1 - nMaxMissingValuesBiomarkers)) & (xControls > 100 * nMaxMissingValuesBiomarkers)) | ...
                ((xPatients > 100 * nMaxMissingValuesBiomarkers) & (xControls < 100 * (1 - nMaxMissingValuesBiomarkers)));
  if sum(pBiomarkers) > 0
    p = find(pBiomarkers > 0);
  else
    p = [];
  end
  % open figure
	hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 
  % plot 
  plot(xPatients, 'o'); hold on; plot(xControls, 'o'); 
  if ~isempty(p)
    plot(p, xPatients(p), 'k.');
    text(p + 50, xPatients(p) + 2, tMetabolites{p});
  end
  line([1, nMetabolites], [100 * nMaxMissingValuesBiomarkers, 100 * nMaxMissingValuesBiomarkers], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
  line([1, nMetabolites], [100 * (1 - nMaxMissingValuesBiomarkers), 100 * (1 - nMaxMissingValuesBiomarkers)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
  box off; xlabel('metabolites'); ylabel('missing values (%)'); xlim([1, nMetabolites]); ylim([0, 100]);
  title(aFile, 'FontWeight', 'normal', 'Interpreter', 'none');
  if ~isempty(p)
    legend({'acromegaly', 'controls', 'biomarkers'});
  else
    legend({'acromegaly', 'controls'});
  end
  % save figure
  aFilename = [aPath, '\\', 'check', '\\', aFile, '_biomarkers.png'];
  print(hFigure, aFilename, '-dpng', '-r300');
  close(hFigure);
  return
end
             
% exclude bad features 
data = data(pBadFeatures == 0, :);

% status
fprintf(1, 'STATUS: %d%% (features) excluded.\n', round(100 * (sum(pBadFeatures == 1) / length(pBadFeatures))));

% sort data
log_data_x = log10(data_x);
x_x = [(-1) * log_data_x(:, labels_x == 1), log_data_x(:, labels_x == 0)];
log_data_y = log10(data);
labels_y = pLabelsClass;
x_y = [(-1) * log_data_y(:, labels_y == 1), log_data_y(:, labels_y == 0)];

% handle data 
x_x = sign(x_x);
y1_x = mean(isnan(x_x(:, labels_x == 1)), 2);
y0_x = mean(isnan(x_x(:, labels_x == 0)), 2);
u1_x = mean(isnan(x_x(:, labels_x == 1)), 1);
u0_x = mean(isnan(x_x(:, labels_x == 0)), 1);
nPatients_x = sum(labels_x == 1);
nControls_x = sum(labels_x == 0);

x_y = sign(x_y);
y1_y = mean(isnan(x_y(:, labels_y == 1)), 2);
y0_y = mean(isnan(x_y(:, labels_y == 0)), 2);
u1_y = mean(isnan(x_y(:, labels_y == 1)), 1);
u0_y = mean(isnan(x_y(:, labels_y == 0)), 1);
nPatients_y = sum(labels_y == 1);
nControls_y = sum(labels_y == 0);

% open figure
hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

% plot original data
subplot(2, 4, [1, 5]); imagesc(x_x, [-4, 4]); colormap('jet'); box off; xlabel('subjects'); ylabel('metabolomics'); xlim([1, nPatients_x + nControls_x]);
title('acromegaly | controls', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
subplot(2, 4, 2); bar([y1_x, y0_x]); box off; ylim([0, 1]); xlabel('metabolomics'); ylabel('proportion of missing values'); 
title('features', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
legend({'acromegaly', 'controls'});
subplot(2, 4, 6); bar(1:nPatients_x, u1_x); hold on; bar((nPatients_x + 1):(nPatients_x + nControls_x), u0_x); box off; ylim([0, 1]); xlabel('subjects'); ylabel('proportion of missing values'); 
title('subjects', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
legend({'acromegaly', 'controls'});

% plot data after exclusion
subplot(2, 4, [3, 7]); imagesc(x_y, [-4, 4]); colormap('jet'); box off; xlabel('subjects'); ylabel('metabolomics'); xlim([1, nPatients_x + nControls_x]); xlim([1, nMetabolites_x]);
title('acromegaly |controls', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
subplot(2, 4, 4); bar([y1_y, y0_y]); box off; ylim([0, 1]); xlabel('metabolomics'); ylabel('proportion of missing values'); 
title('features', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
legend({'acromegaly', 'controls'});
subplot(2, 4, 8); bar(1:nPatients_y, u1_y); hold on; bar((nPatients_y + 1):(nPatients_y + nControls_y), u0_y); box off; ylim([0, 1]); xlabel('subjects'); ylabel('proportion of missing values'); 
title('subjects', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
legend({'acromegaly', 'controls'});

% save figure
aFile = sprintf('%s_MV%d%d', aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
aFilename = [aPath, '\\', 'check', '\\', aFile, '_exclusion.png'];
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

% save
bSaveLucyData = 0;
if bSaveLucyData == 1
  x = [nanmean(log_data_y(:, labels_y == 1), 2), nanmean(log_data_y(:, labels_y == 0), 2)];
  % convert to cells
  nFeatures = size(data, 1);
  tCell = cell(nFeatures, 3); 
  % init
  tCell(1:end, 1) = names(pBadFeatures == 0);
  tCell(1:end, 2:3) = num2cell(x);
  % write table
  tTable = cell2table(tCell, 'VariableNames', {'Class', 'Acromegaly', 'Controls'});
  aFilename = [aPath, '\\', 'check', '\\', aFile, '_mu.csv'];
  writetable(tTable, aFilename);
end

end % end

%-------------------------------------------------------------------------------