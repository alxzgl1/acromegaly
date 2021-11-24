%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_import()

clc;

aPath = 'd:\\data\\acromegaly\\import';
aFile = 'Lipids_P'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

% parameters
nMaxMissingValuesBySubjects = 0.30;
nMaxMissingValuesByFeatures = 0.50; % keep above 90% separately as potential biomarker

% load data and header
aFilename = [aPath, '\\', aFile, '_data.mat'];
load(aFilename, 'data');
aFilename = [aPath, '\\', aFile, '_names.mat'];
load(aFilename, 'names');
aFilename = [aPath, '\\', aFile, '_labels.mat'];
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

% process data | exclude missing values by subject and by feature
% * by subject: exclude subject if 30% of data are missing
% * by feature: exclude feature if 50% for both groups are missing
%               keep separately if 90% missing in one of the group
         
% count missing values by subjects       
nFeatures = size(data, 1);
pMissingValuesBySubjects = sum(isnan(data)) / nFeatures;
pBadSubjects = pMissingValuesBySubjects > nMaxMissingValuesBySubjects;

% exclude bad subjects and labels
data = data(:, pBadSubjects == 0);
tLabelsID = tLabelsID(pBadSubjects == 0);
tLabelsClass = tLabelsClass(pBadSubjects == 0);
pLabelsClass = pLabelsClass(pBadSubjects == 0);
tLabelsSex = tLabelsSex(pBadSubjects == 0);
tLabelsAge = tLabelsAge(pBadSubjects == 0);
tLabelsBMI = tLabelsBMI(pBadSubjects == 0);
 
% count missing values by features
nPatients = sum(pLabelsClass == 1);
nControls = sum(pLabelsClass == 0);
pMissingValuesByFeatures = [sum(isnan(data(:, pLabelsClass == 1)), 2) / nPatients, ...
                            sum(isnan(data(:, pLabelsClass == 0)), 2) / nControls];

pBadFeatures = pMissingValuesByFeatures(:, 1) > nMaxMissingValuesByFeatures | ...
               pMissingValuesByFeatures(:, 2) > nMaxMissingValuesByFeatures;
             


% convert to cells
nSubjects = size(data, 2);
nFeatures = size(data, 1);
tCell = cell(nFeatures + 1, nSubjects + 1); 

% fill info
tCell(1, 1:end) = ['Class', tLabelsClass];
tCell(2:end, 1) = names;

% fill data
for iFeature = 1:nFeatures
  for iSubject = 1:nSubjects
    tCell{iFeature + 1, iSubject + 1} = data(iFeature, iSubject);
  end
end

tTable = cell2table(tCell, 'VariableNames', ['ID', tLabelsID]);

% write table
writetable(tTable, [aFile, '.csv']);

% % open figure
% hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.5, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 
% 
% % plot
% subplot(2, 2, [1, 3]); imagesc(x, [-4, 4]); colormap('jet'); box off; xlabel('subjects'); ylabel('metabolomics'); 
% title('Missing values | acromegaly (cyan) and controls (yellow)', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
% subplot(2, 2, 2); bar([y1, y0]); box off; ylim([0, 1]); xlabel('metabolomics'); ylabel('proportion of missing values'); 
% title('Missing values | acromegaly (blue) and controls (orange)', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
% subplot(2, 2, 4); bar(find(labels == 1), u1); hold on; bar(find(labels == 0), u0); box off; ylim([0, 1]); xlabel('subjects'); ylabel('proportion of missing values'); 
% title('Missing values | acromegaly (blue) and controls (orange)', 'FontWeight', 'normal', 'FontSize', 8); set(gca, 'FontSize', 8);
% legend({'acromegaly', 'controls'});
% 
% % save figure
% aFilename = [aPath, \\', aFile, '.png'];
% print(hFigure, aFilename, '-dpng', '-r300');
% close(hFigure);

end % end

%-------------------------------------------------------------------------------