% Function that calculates the current cost statistic of the given data and prototypes/omega
% @param this GMLVQ.Run
% @param dataPair GMLVQ.DataPair
% @param prototypes
% @param omegaMatrix
% @param customMu double optional
% @out costFunction double
% @out crispOut 1xN double vector
% @out margins 1xN double vector
% @out score 1xN double vector
% (with N = nFeatureVectors in the dataPair input)
function [costFunction, crispOut, margins, score] = cost(this, dataPair, prototypes, omegaMatrix, customMu)

if nargin < 5
    mu = this.gmlvq.params.mu;
else
    mu = customMu;
end

% Initialization
costFunction = 0;
margins = zeros(1, dataPair.nFeatureVectors);
score = zeros(1, dataPair.nFeatureVectors);
crispOut = zeros(1, dataPair.nFeatureVectors);

omegaMatrix = omegaMatrix / norm(omegaMatrix,'fro');

% Loop through all examples to find the cost of each
for i=1:dataPair.nFeatureVectors
    currentFV = dataPair.featureVectors(i,:);
    currentLabel = dataPair.labels(i);
    dist = nan(this.nPrototypes, 1);
    
    for jk = 1:this.nPrototypes
        dist(jk) = norm(omegaMatrix * (currentFV-prototypes(jk,:))')^2;
    end
    
    % Find the winning prototypes for this example
    correct = find(this.gmlvq.plbl == currentLabel);
    incorrect = find(this.gmlvq.plbl ~= currentLabel);
    [dJJ,JJJ] = min(dist(correct));
    [dKK,KKK] = min(dist(incorrect));
    JJ = correct(JJJ); KK = incorrect(KKK);
    
    margins(i) = (dJJ-dKK)/(dJJ+dKK);
    costFunction = costFunction + margins(i) / dataPair.nFeatureVectors;
    
    % Non-normalized difference of distances
    if currentLabel == this.gmlvq.params.rocClass
        score(i) = dKK - dJJ;
    else
        score(i) = dJJ - dKK;
    end
    
    % Class label according to the nearest prototype
    crispOut(i) = this.gmlvq.plbl(JJ) * (dJJ <= dKK) + this.gmlvq.plbl(KK) * (dJJ > dKK);
end

% Penalty term
if mu > 0
    costFunction = costFunction - mu / ...
        2 * log(det(omegaMatrix*omegaMatrix')) / dataPair.nFeatureVectors;
end

end