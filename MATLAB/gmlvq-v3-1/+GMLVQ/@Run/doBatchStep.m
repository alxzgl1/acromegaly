% Function that does the real step in the algorithm
% @param this GMLVQ.Run
% @param prototypes
% @param omegaMatrix
% @out prototypes
% @out omegaMatrix
function [prototypes, omegaMatrix] = doBatchStep(this, prototypes, omegaMatrix)

lambda = omegaMatrix'*omegaMatrix;

% Initialize change in prot, omega
changePrototypes = 0*prototypes;
changeOmega = 0*omegaMatrix;

% Loop through all training samples
for i = 1:this.trainingData.nFeatureVectors
    currentFV = this.trainingData.featureVectors(i,:);
    currentLabel = this.trainingData.labels(i);
    
    % Calculate squared distances to all prototypes
    dist = nan(this.nPrototypes, 1);
    for j = 1:this.nPrototypes
        dist(j) = (norm(omegaMatrix * (currentFV - prototypes(j,:))'))^2;
    end
    
    % Find the winning prototypes
    correct = find(this.gmlvq.plbl == currentLabel);
    incorrect = find(this.gmlvq.plbl ~= currentLabel);
    
    [dJ, JJ] = min(dist(correct));
    [dK,KK] = min(dist(incorrect));
    
    % Get winning prototypes and indices
    jwin = correct(JJ);
    kwin = incorrect(KK);
    wJ = prototypes(jwin, :);
    wK = prototypes(kwin, :);
    
    % GMLVQ prototype update for one example fvi
    DJ = (currentFV - wJ)'; % displacement vector
    DK = (currentFV - wK)';
    normalizationFactor = (dJ + dK)^2;
    dwJ = -(dK/normalizationFactor)*lambda*DJ;
    dwK = (dJ/normalizationFactor)*lambda*DK;
    
    % Matrix update, single matrix omegaMatrix for one example
    f1 = (dK/normalizationFactor)*(omegaMatrix*DJ)*DJ';
    f2 = (-dJ/normalizationFactor)*(omegaMatrix*DK)*DK';
    
    % Negative gradient added up over examples
    changePrototypes(jwin,:) = changePrototypes(jwin,:) - dwJ';
    changePrototypes(kwin,:) = changePrototypes(kwin,:) - dwK';
    changeOmega = changeOmega - (f1 + f2);
end

% Singularity control: Add derivative of penalty term times mu
if this.gmlvq.params.mu > 0
    changeOmega = changeOmega + this.gmlvq.params.mu * pinv(omegaMatrix');
end

% Compute normalized gradient updates
% Separate normalization for prototypes and the matrix
% Computation of actual changes, diagonal matrix imposed here if necessary
if this.gmlvq.params.mode == GMLVQ.Mode.GRLVQ
    changeOmega = diag(diag(changeOmega));
end

% Final, normalized gradient updates after 1 loop through examples
prototypes = prototypes + this.gmlvq.params.etap * changePrototypes / norm(changePrototypes,'fro');
omegaMatrix = omegaMatrix + this.gmlvq.params.etam * changeOmega / norm(changeOmega,'fro');

if this.gmlvq.params.mode == GMLVQ.Mode.GRLVQ
    omegaMatrix = diag(diag(omegaMatrix));
end

% Null-space correction using Moore-Penrose pseudo-inverse
if this.gmlvq.params.mode == GMLVQ.Mode.GMLVQNS
    xvec = [this.trainingData.featureVectors; prototypes];
    omegaMatrix = (omegaMatrix * xvec') * pinv(xvec');
end

if this.gmlvq.params.mode == GMLVQ.Mode.GLVQ
    omegaMatrix = eye(this.nDimensions);
end

% Normalization of omega, corresponds to Trace(lambda) = 1
omegaMatrix = omegaMatrix / norm(omegaMatrix,'fro');

end