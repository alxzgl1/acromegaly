% Initialization of the prototypes close to class conditional means, small random displacements
% applied to break ties.
% @param this GMLVQ.Run
% @out prototypes
% @out omegaMatrix
function [prototypes, omegaMatrix] = setInitialPrototypes(this)

% RJV Repeat warnings and checks here
this.gmlvq.params.sanityCheck();

prototypes = zeros(this.nPrototypes, this.nDimensions);

% Compute class-conditional means
for i = 1:this.nPrototypes
    prototypes(i,:) = ...
        mean(this.trainingData.featureVectors(this.trainingData.labels == this.gmlvq.plbl(i),:),1);
end

% Random displacement
prototypes = prototypes .* (0.99 + this.gmlvq.params.randomization*rand(size(prototypes)));

% Run classwise k-means
if this.gmlvq.params.useKMeans
    for i = 1:this.nClasses
        stream = RandStream('mlfg6331_64');
        options = statset('UseParallel', 1, 'UseSubstreams', 1, 'Streams', stream);
        noPrototypes = sum(this.gmlvq.plbl == i);

        % RJV Editor's Note: kmeans will automatically use Parallel Toolbox
        % if installed.
        [~, prototypes(this.gmlvq.plbl == i,:)] = kmeans(this.trainingData.featureVectors(this.trainingData.labels == i,:), noPrototypes,...
            'Options', options, 'Replicates', 5);
    end
end

% Matrix initialization, identity or random
omegaMatrix = eye(this.nDimensions); % this works if rndinit param == 0

% Fix for GLVQ (mode 3)
if this.gmlvq.params.mode ~= GMLVQ.Mode.GLVQ && this.gmlvq.params.rndinit
    omegaMatrix = rand(this.nDimensions) - 0.5;
    omegaMatrix = omegaMatrix' * omegaMatrix;
end

% For GRLVQ we restrict to a diagonal matrix
if this.gmlvq.params.mode == GMLVQ.Mode.GRLVQ
    omegaMatrix = diag(diag(omegaMatrix));
end

% Normalization, to make trace(lambda) = 1
omegaMatrix = omegaMatrix/sqrt(sum(sum(omegaMatrix.^2)));

end