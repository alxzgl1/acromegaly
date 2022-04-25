function checkConsistency(gmlvq)

% General parameter analysis
if gmlvq.data.nFeatureVectors ~= length(gmlvq.data.labels)
    error('Number of training labels differs from number of samples!');
end

if gmlvq.nClasses > 2
    warning off backtrace
    warning(['Multi-class problem. ROC analysis is for class ', num2str(gmlvq.params.rocClass), ' (neg.) vs. all others (pos.)']);
    warning on backtrace
end

if length(unique(gmlvq.plbl)) ~= length(unique(gmlvq.data.labels))
    error('Number of prototype classes must equal number of classes!');
end

if sum(unique(gmlvq.plbl) ~= unique(gmlvq.data.labels)) > 0
    error('Prototype labels are inconsistent with data, please rename or reorder');
end

% Standard deviation analysis
sd = std(gmlvq.data.featureVectors, 0, 1);
disp(' ');
disp(['Minimum standard deviation of features: ', num2str(min(sd))]);
if min(sd) < 1.e-10
    error('At least one feature shows (close to) zero variance');
end

if gmlvq.params.ncop >= gmlvq.totalsteps
    error('Number of gradient steps must be larger than ncop!');
end

end