function h = plotconf(this, type, stepID, suffix)
%PLOTCONF Plot a confusion matrix
%   Use this to plot a training or validation matrix from a Run or
%   averageRun
if nargin < 2
    type = 'training';
end

if ~any(strcmp({'training','validation'},type))
   error('type must be either ''training'' or ''validation'''); 
end

if nargin > 3 
    suffix = 'Std';
end

if nargin < 4 || isempty(suffix), suffix = ''; end

dataprop = [type 'Perf' suffix];

if (nargin < 3) || (stepID == -1)
    stepID = length(this.(dataprop));
end

if stepID == 0 || isempty(this.(dataprop))
   error('Plot %s not available in this run', dataprop);
end

confmat = this.(dataprop)(stepID).confusionMatrix;

% figure;
h = GMLVQ.Util.confplot_rjv(confmat, sprintf('%s(%i)', dataprop, stepID));

end