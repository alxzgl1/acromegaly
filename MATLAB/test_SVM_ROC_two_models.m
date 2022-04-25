%-------------------------------------------------------------------------------
% Function 
% Reference: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/65629/versions/1/previews/plot_roc_curve_example_comparealgos.m/index.html
% https://uk.mathworks.com/help/stats/perfcurve.html
%-------------------------------------------------------------------------------
function test_SVM_ROC_two_models()

clc;
rng(1);  % For reproducibility
n = 100; % Number of points per quadrant

r1 = sqrt(rand(2*n,1));                     % Random radii
t1 = [pi/2*rand(n,1); (pi/2*rand(n,1)+pi)]; % Random angles for Q1 and Q3
X1 = [r1.*cos(t1) r1.*sin(t1)];             % Polar-to-Cartesian conversion

r2 = sqrt(rand(2*n,1));
t2 = [pi/2*rand(n,1)+pi/2; (pi/2*rand(n,1)-pi/2)]; % Random angles for Q2 and Q4
X2 = [r2.*cos(t2) r2.*sin(t2)];

pred = [X1; X2];
resp = ones(4*n, 1);
resp(2*n + 1:end) = -1; % Labels

SVMModel1 = fitcsvm(pred, resp, 'KernelFunction', 'mysigmoid_A', 'Standardize', true);
SVMModel1 = fitPosterior(SVMModel1);
[~, scores1] = resubPredict(SVMModel1);

SVMModel2 = fitcsvm(pred, resp, 'KernelFunction', 'mysigmoid_B', 'Standardize', true);
SVMModel2 = fitPosterior(SVMModel2);
[~, scores2] = resubPredict(SVMModel2);

[x1, y1, ~, auc1] = perfcurve(resp, scores1(:, 2), 1);
[x2, y2, ~, auc2] = perfcurve(resp, scores2(:, 2), 1);

plot(x1,y1)
hold on
plot(x2,y2)
hold off
legend('gamma = 1','gamma = 0.5','Location','SE');
xlabel('False positive rate'); ylabel('True positive rate');
title('ROC for classification by SVM');

end % end

%-------------------------------------------------------------------------------