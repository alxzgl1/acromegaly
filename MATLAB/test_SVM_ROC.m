%-------------------------------------------------------------------------------
% Function 
% Reference: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/65629/versions/1/previews/plot_roc_curve_example_comparealgos.m/index.html
%-------------------------------------------------------------------------------
function test_SVM_ROC()

clc;

aModel = 'SVM'; % 'SVM', 'GLM'

load fisheriris

pred = meas(51:end,1:2);

resp = (1:100)' > 50;  % Versicolor = 0, virginica = 1 | binary response variable

% model
if strcmp(aModel, 'GLM') % fit a logistic regression model
  mdl = fitglm(pred, resp, 'Distribution', 'binomial', 'Link', 'logit');
  % Compute the ROC curve. Use the probability estimates from the logistic regression model as scores.
  scores = mdl.Fitted.Probability;
  [X, Y, T, AUC] = perfcurve(species(51:end, :), scores, 'virginica');
elseif strcmp(aModel, 'SVM')
  mdl = fitcsvm(pred, resp, 'KernelFunction', 'linear', 'KernelScale', 'auto');
  % compute the posterior probabilities (scores)
  mdl = fitPosterior(mdl);
  [~, scores] = resubPredict(mdl);
  % compute the standard ROC curve using the scores from the SVM model.
  [X, Y, T, AUC] = perfcurve(species(51:end, :), scores(:, mdl.ClassNames), 'virginica');
end

% display the area under the curve
AUC

% plot the ROC curve
plot(X, Y);
xlabel('False positive rate');
ylabel('True positive rate');
title(sprintf('ROC for Classification by %s', aModel));

end % end

%-------------------------------------------------------------------------------