%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function [w, accuracy] = support_average_models(pAccuracies, pCoeffs, aModelAveraging)

pWeights = pAccuracies;
nModels = length(pWeights);

% maximum / weighted
if strcmp(aModelAveraging, 'mean') 
  w = mean(pCoeffs, 2);
  accuracy = mean(pAccuracies);
elseif strcmp(aModelAveraging(1:4), 'poly')
  % scale weights
  degree = str2double(aModelAveraging(5));
  pWeights = nModels * ((pWeights .^ degree) / sum(pWeights .^ degree));
  % weighted coeffs
  nCoeffs = size(pCoeffs, 1);
  w = mean(pCoeffs .* repmat(pWeights', nCoeffs, 1), 2);
  % accuracy
  accuracy = mean(pAccuracies .* pWeights);
end

end % end

%-------------------------------------------------------------------------------