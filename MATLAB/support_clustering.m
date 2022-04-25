%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function [pA, pH, nThreshold, pI] = support_clustering(pA, nDepth, bDepth)

% clustering
% opt.type_sim = 'average', 'single', 'complete', 'average', 'ward' (default)
opt = [];
pH = support_hierarchical_clustering(pA, opt);
pI = [];

% unbundle hierarchy
nThreshold = -1;
if nDepth ~= -1
  nItems = pH(1, 4) - 1;
  pLevel = pH(:, 1);
  
  if bDepth == 1
    nThreshold = (pLevel(1) - pLevel(end)) * (1 - nDepth) + pLevel(end);
  else
    nClusters = nDepth;
    nThreshold = pLevel(end - (nClusters - 1));
  end
  
  pLevel(pLevel < nThreshold) = [];
  nLevels = length(pLevel);
  pX = pH(1:nLevels, 2);
  pY = pH(1:nLevels, 3);
  pZ = pH(1:nLevels, 4);
  pI = (1:nItems); pI = pI(:);
  for n = 1:nLevels
    pXIndex = find(pI == pX(n));
    pYIndex = find(pI == pY(n));
    if ~isempty(pXIndex), pI(pXIndex) = pZ(n); end
    if ~isempty(pYIndex), pI(pYIndex) = pZ(n); end
    if pX(n) < (nItems + 1), pI(pX(n)) = pZ(n); end
    if pY(n) < (nItems + 1), pI(pY(n)) = pZ(n); end
  end
  % get unique clusters
  pUnique = unique(pI);
  for n = 1:length(pUnique)
    pI(pI == pUnique(n)) = n;
  end
  for n = 1:length(pI)
    for m = 1:length(pI)
      pA(n, m) = 0;
      if pI(n) == pI(m), pA(n, m) = pI(n); end
    end
  end
end

end %% end

%-------------------------------------------------------------------------------