%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_data_import()

clc;

aPath = 'd:\\data\\acromegaly\\import';
aFile = 'Lipids_P'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

% note: biomarker 'HILIC_P' is detected ('M235T108')

% parameters
nMaxMissingValuesBySubjects = 0.30;
nMaxMissingValuesByFeatures = 0.50; 
nMaxMissingValuesBiomarkers = 0.80; % keep above 80% separately as potential biomarker

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
         
% count missing values by subjects       
nFeatures = size(data, 1);
pMissingValuesBySubjects = sum(isnan(data)) / nFeatures;
pBadSubjects = pMissingValuesBySubjects > nMaxMissingValuesBySubjects;

nMaxPatients = sum(pLabelsClass == 1);
nMaxControls = sum(pLabelsClass == 0);

% exclude bad subjects 
data = data(:, pBadSubjects == 0);
tLabelsID = tLabelsID(pBadSubjects == 0);
tLabelsClass = tLabelsClass(pBadSubjects == 0);
pLabelsClass = pLabelsClass(pBadSubjects == 0);
tLabelsSex = tLabelsSex(pBadSubjects == 0);
tLabelsAge = tLabelsAge(pBadSubjects == 0);
tLabelsBMI = tLabelsBMI(pBadSubjects == 0);

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

% exclude bad features 
data = data(pBadFeatures == 0, :);
names = names(pBadFeatures == 0);

% status
fprintf(1, 'STATUS: %d%% (features) excluded.\n', round(100 * (sum(pBadFeatures == 1) / length(pBadFeatures))));
 
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
aFilename = sprintf('%s\\%s_MV%d%d.csv', aPath, aFile, round(100 * nMaxMissingValuesBySubjects), round(100 * nMaxMissingValuesByFeatures));
writetable(tTable, aFilename);

end % end

%-------------------------------------------------------------------------------