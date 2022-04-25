%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function test_display_metabolites()

clc;

aFile = 'HILIC_P'; % 'HILIC_N', HILIC_P', 'Lipids_N', 'Lipids_P'

nSplit = 0.70; 

pEVs = [0.60, 0.70, 0.75, 0.80, 0.90, 0.95, 0.99];
nEVs = length(pEVs); 

nMaxMetabolites = 200;

pMetabolites = zeros(nMaxMetabolites, nEVs);

% load
for iEV = 1:nEVs
  nEV = pEVs(iEV);
  % load metabolites
  aFilename = sprintf('d:\\data\\acromegaly\\_analysis\\PCA_SVM\\metabolites\\%s_%d%d_EV_%d.mat', aFile, 100 * nSplit, round(100 * (1 - nSplit)), round(100 * nEV));
  load(aFilename, 'pMetaboliteIndices', 'names');
  % init
  pMetabolites(:, iEV) = pMetaboliteIndices(end:-1:(end - nMaxMetabolites + 1));
end

X = pMetabolites(:);
pU = unique(X);

nU = length(pU);

pRanks = zeros(nU, nEVs);
pRanksNonzero = zeros(nU, nEVs);

for iU = 1:nU
  u = pU(iU);
  for iEV = 1:nEVs
    x = pMetabolites(:, iEV);
    i = find(x == u);
    if ~isempty(i)
      pRanks(iU, iEV) = i;
      pRanksNonzero(iU, iEV) = pRanksNonzero(iU, iEV) + 1;
    else
     pRanks(iU, iEV) = nMaxMetabolites;
    end
  end
end

pRanksNonzeroMu = mean(pRanksNonzero, 2);

pRanksMu = mean(pRanks, 2);
[~, iRanks] = sort(pRanksMu);

pU_Ranks = pU(iRanks);

for i = 1:nMaxMetabolites
  fprintf(1, '%d\t%1.2f\t%1.2f\t%s\n', pU_Ranks(i), pRanksMu(iRanks(i)), pRanksNonzeroMu(iRanks(i)), names{pU_Ranks(i)});
end

% display
% nMaxMetabolites = 100;
% for i = 1:nMaxMetabolites
%   fprintf(1, '%d\n', pMetaboliteIndices(end - i + 1));
% end

end % end

%-------------------------------------------------------------------------------