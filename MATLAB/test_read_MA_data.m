%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function test_read_MA_data()

clc;

aPath = 'd:\\data\\acromegaly\\import\\MA\\univariate';
aFile = 'Lipids_N'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

% parameters
nMV_BySubjects = 0.30;
nMV_ByFeatures = 0.25; % controls should be equal to patients | 0.50, 0.25

% load data
aFilename = sprintf('%s\\%s_S%dF%d_NA.csv', aPath, aFile, round(100 * nMV_BySubjects), round(100 * nMV_ByFeatures));
T = readcell(aFilename);

% parse data
data = cell2mat(T(3:end, 2:end));
labels = contains(T(2, 2:end), 'Acromegaly');
names = T(3:end, 1);

end % end

%-------------------------------------------------------------------------------