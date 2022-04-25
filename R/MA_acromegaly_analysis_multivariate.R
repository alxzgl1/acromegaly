# load package
library("MetaboAnalystR")

# set directory
setwd("D:/code/__R/Seedcorn") # L15 laptop

# filename
aFilename = "D:/data/acromegaly/import/HILIC_N_S30F25_NA.csv" # HILIC_N, HILIC_P, Lipids_N, Lipids_P

# init
mSet = InitDataObjects("pktable", "stat", FALSE)
mSet = Read.TextData(mSet, aFilename, "colu", "disc") 

# impute values
mSet = ImputeMissingVar(mSet, method="knn_var")

# WARNING: SanityCheckData
#  SanityCheckData is used for data processing, and performs a basic sanity check of the
#  uploaded content, ensuring that the data is suitable for further analysis. The function will
#  return a message if the data has successfully passed the check and is deemed suitable for
#  further analysis. If it fails, the function will return a 0. The function will perform the
#  check directly onto the mSet$dataSet object, and must be performed immediately after
#  reading in data. The sanity check function evaluates the accuracy of sample and class
#  labels, data structure, deals with non-numeric values, removes columns that are constant
#  across all samples (variance = 0), and by default replaces missing values with half of the
#  original minimal positive value in your dataset.

# prenorm
mSet = PreparePrenormData(mSet)

# normalisation + log10 + pareto scaling 
mSet = Normalization(mSet, "SumNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
# mSet = PlotNormSummary(mSet, "multivariate_norm_4_", "png", 72, width=NA)
# mSet = PlotSampleNormSummary(mSet, "multivariate_snorm_4_", "png", 72, width=NA)

# sanity check
mSet = SanityCheckData(mSet)

# Correlation Heatmaps
bCORMAP = FALSE
if (bCORMAP == TRUE) {
  mSet = PlotCorrHeatMap(mSet, "corr_0_", "png", 72, width=NA, "col", "pearson", "bwm", "overview", F, F, "0")
}

# PCA
bPCA = TRUE
if (bPCA == TRUE) {
  mSet = PCA.Anal(mSet)
  mSet = PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
  mSet = PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
  mSet = PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1, 2, 0.95, 1, 0)
  mSet = PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1, 2);
  mSet = PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1, 2)
  mSet = PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1, 2, 3)
}

# Partial Least Squares Discriminant Analysis (PLS-DA)
bPLS_DA = TRUE
if (bPLS_DA == TRUE) {
  mSet = PLSR.Anal(mSet, reg=TRUE)
  mSet = PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5)
  mSet = PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
  mSet = PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
  mSet = PlotPLSLoading(mSet, "pls_loading_0_", "png", 72, width=NA, 1, 2);
  mSet = PlotPLS3DLoading(mSet, "pls_loading3d_0_", "json", 1,2,3)
  mSet = PLSDA.CV(mSet, "T", 5, "Q2")
  mSet = PlotPLS.Classification(mSet, "pls_cv_0_", "png", 72, width=NA)
  mSet = PlotPLS.Imp(mSet, "pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)
}

# Sparse PLS-DA
bSPLS_DA = TRUE
if (bSPLS_DA == TRUE) {
  mSet = SPLSR.Anal(mSet, 5, 10, "same", "Mfold")
  mSet = PlotSPLSPairSummary(mSet, "spls_pair_0_", "png", 72, width=NA, 5)
  mSet = PlotSPLS2DScore(mSet, "spls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
  mSet = PlotSPLS3DScoreImg(mSet, "spls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
  mSet = PlotSPLSLoading(mSet, "spls_loading_0_", "png", 72, width=NA, 1,"overview");
  mSet = PlotSPLSDA.Classification(mSet, "spls_cv_0_", "png", 72, width=NA)
  mSet = PlotSPLS3DLoading(mSet, "spls_loading3d_0_", "json", 1,2,3)
}

# Orthogonal PLS-DA 
bOPLS_DA = TRUE
if (bOPLS_DA == TRUE) {
  mSet = OPLSR.Anal(mSet, reg=TRUE)
  mSet = PlotOPLS2DScore(mSet, "opls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
  mSet = PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
  mSet = PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE)
  mSet = PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)
}

# Random Forest
bRF = TRUE
if (bRF == TRUE) {
  mSet = RF.Anal(mSet, 500,7,1)
  mSet = PlotRF.Classify(mSet, "rf_cls_0_", "png", 72, width=NA)
  mSet = PlotRF.VIP(mSet, "rf_imp_0_", "png", 72, width=NA)
  mSet = PlotRF.Outlier(mSet, "rf_outlier_0_", "png", 72, width=NA)
}

# SVM
bSVM = TRUE
if (bSVM == TRUE) {
  mSet = RSVM.Anal(mSet, 10)
  mSet = PlotRSVM.Classification(mSet, "svm_cls_0_", "png", 72, width=NA)
  mSet = PlotRSVM.Cmpd(mSet, "svm_imp_0_", "png", 72, width=NA)
}
