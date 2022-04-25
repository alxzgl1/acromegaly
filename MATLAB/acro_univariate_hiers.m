%-------------------------------------------------------------------------------
% Function 
%-------------------------------------------------------------------------------
function acro_univariate_hiers()

clc;

aFile = 'Lipids_N'; % 'HILIC_N', 'HILIC_P', 'Lipids_N', 'Lipids_P'

aAnalysisType = 'univariate'; % 1/5 min imputation, norm by sum, log10, no-scaling

nMinPValue = 0.001; 
nClusters = 4; % number of clusters | adjust after visual check

% parameters
nMV_BySubjects = 0.30;
nMV_ByFeatures = 0.25; % 0.25 (default)

% path
aPath = 'd:\\data\\acromegaly';

% open figure
hFigure = figure; set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 1.5, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

% load data
aExclusion = sprintf('S%dF%d', round(100 * nMV_BySubjects), round(100 * nMV_ByFeatures));
aFilename = sprintf('%s\\import\\MA\\%s\\%s_%s_NA.csv', aPath, aAnalysisType, aFile, aExclusion);
T = readcell(aFilename);
% parse data
data = cell2mat(T(3:end, 2:end)); 
pLabelsClass = contains(T(2, 2:end), 'Acromegaly');
names = T(3:end, 1);

% init 
x = data;
tMetabolites = names;

% groups
X0 = x(:, pLabelsClass == 0);
X1 = x(:, pLabelsClass == 1);

% loop features
nFeatures = size(x, 1);
pPValues = zeros(nFeatures, 1);
pPValuesRaw = zeros(nFeatures, 1);
for iFeature = 1:nFeatures
  x0 = X0(iFeature, :)';
  x1 = X1(iFeature, :)';
  bParametric = 0;
  if bParametric == 1
    [~, pval] = ttest2(x0, x1, 'Vartype', 'unequal');
  else
    pval = ranksum(x0, x1);
  end
  pPValues(iFeature) = -log10(pval);
  pPValuesRaw(iFeature) = pval;
end

% FDR Benjamini-Hochberg
bhFDR = fdr_bh(pPValuesRaw, 0.05);
fprintf(1, 'Significant p-values after BH-FDR: %d\n', sum(bhFDR));
% FDR Storey
[sFDR, sQ] = mafdr(pPValuesRaw);
sFDR = sFDR < 0.05;
sQ = sQ < 0.05;
fprintf(1, 'Significant p-values after S-FDR: %d\n', sum(sFDR));
fprintf(1, 'Significant p-values after S-Q: %d\n', sum(sQ));
fprintf(1, 'BH-FDR & S-FDR: %d\n', sum(bhFDR & sFDR));
fprintf(1, 'BH-FDR & S-Q: %d\n', sum(bhFDR & sQ));

[px, ix] = sort(pPValues, 'descend');
nThresholds = find(px > -log10(nMinPValue));
if isempty(nThresholds)
  return
end
tMetabolites = tMetabolites(ix);
pPValues = pPValues(ix);
bhFDR = bhFDR(ix);
sFDR = sFDR(ix);
sQ = sQ(ix);
nThresholds = nThresholds(end);
ix = ix(1:nThresholds);

y = x(ix, :);

tMetabolites = tMetabolites(1:nThresholds);
pPValues = pPValues(1:nThresholds);
bhFDR = bhFDR(1:nThresholds);
sFDR = sFDR(1:nThresholds);
sQ = sQ(1:nThresholds);

% correlation
CM = corr(y', 'type', 'Spearman'); % Peter suggested abs()

% clutering
[pA, pH, nDepth, pI] = support_clustering(CM, nClusters, 0);

% plot
subplot(2, 3, 1); imagesc(CM, [-1, 1]); set(gca, 'YDir', 'normal'); colorbar;
xlabel('metabolites'); ylabel('metabolites'); title(sprintf('CM | %s', aFile), 'FontWeight', 'normal', 'Interpreter', 'none'); box off;

subplot(2, 3, 2); barh(pPValues); box off; ylim([0.5, nThresholds + 0.5]); hold on;
line([-log10(nMinPValue), -log10(nMinPValue)], [0.5, nThresholds + 0.5], 'Color', 'g', 'LineWidth', 2);
if sum(bhFDR) > 0
  nFDR_PValue = pPValues(find(bhFDR == 0, 1) - 1);
  line([nFDR_PValue, nFDR_PValue], [0.5, nThresholds + 0.5], 'Color', 'r', 'LineWidth', 2);
end
if sum(sFDR) > 0
  nFDR_PValue = pPValues(find(sFDR == 0, 1) - 1);
  line([nFDR_PValue, nFDR_PValue], [0.5, nThresholds + 0.5], 'Color', 'm', 'LineWidth', 2);
end
if sum(sQ) > 0
  nFDR_PValue = pPValues(find(sQ == 0, 1) - 1);
  line([nFDR_PValue, nFDR_PValue], [0.5, nThresholds + 0.5], 'Color', 'c', 'LineWidth', 2, 'LineStyle', '--');
end
ylabel('metabolites'); xlabel('-log10(pval)');
title(sprintf('%s | BH-FDR (r) | S-FDR/S-Q (m/c)', aFile), 'FontWeight', 'normal', 'Interpreter', 'none');

subplot(2, 3, 3); support_dendrogram(pH); title(sprintf('dendrogram | depth=%1.3f', nDepth), 'FontWeight', 'normal');
xlabel('metabolites'); 

subplot(2, 3, 6); imagesc(pA); set(gca, 'YDir', 'normal'); title('clusters', 'FontWeight', 'normal'); box off;
xlabel('metabolites'); ylabel('metabolites'); 

subplot(2, 3, 5); barh(pI); title('node labels', 'FontWeight', 'normal'); box off;
ylim([0.5, nThresholds + 0.5]);
ylabel('metabolites'); xlabel('clusters id');

% save figure
aFilename = [aPath, '\\', '_analysis', '\\', [aAnalysisType, '_hiers'], '\\', [aFile, '_', aExclusion], '.png'];
print(hFigure, aFilename, '-dpng', '-r300');
close(hFigure);

% sort pI
[~, h] = sort(pI);

% save names
hFile = fopen([aFilename(1:(end - 4)), '.txt'], 'wt');
for i = 1:length(pI)
  fprintf(hFile, '%d\t%s\t%1.3f\t%d\t%d\t%d\n', pI(h(i)), tMetabolites{h(i)}, pPValues(h(i)), bhFDR(h(i)), sFDR(h(i)), sQ(h(i)));
end
fclose(hFile);

end % end

%-------------------------------------------------------------------------------