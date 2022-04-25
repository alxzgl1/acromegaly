%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function test_display_w()

clc;

aFile = 'Lipids_N'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

nSplit = 0.50; 

pEVs = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95];
nEVs = length(pEVs);

aModelAveraging = 'poly2';

nMaxSubjects = 100;

X = zeros(nMaxSubjects, nEVs);
pAccuracies = zeros();

% load
for iEV = 1:nEVs
  nEV = pEVs(iEV);
  % load w
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\weights\\%s_%d%d_EV_%d_%s.mat', aFile, 100 * nSplit, round(100 * (1 - nSplit)), round(100 * nEV), aModelAveraging);
  load(aFilename, 'w', 'x', 'accuracy', 'subjects', 'labels'); % 'K'

  % init
  nSubjects = length(x);
  X(1:nSubjects, iEV) = x;

  pAccuracies(iEV) = accuracy;
end

X = X(1:nSubjects, :);

for iEV = 1:nEVs
  for iSubject = 1:nSubjects
    fprintf(1, '%1.4f\t', X(iSubject, iEV));
  end
  fprintf(1, '\n');
end

for iSubject = 1:nSubjects
  fprintf(1, '%s\t', subjects{iSubject});
end
fprintf(1, '\n');

end % end

%-------------------------------------------------------------------------------