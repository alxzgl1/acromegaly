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
mSet = ImputeMissingVar(mSet, method="min")

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

# pernorm
mSet = PreparePrenormData(mSet)

# normalisation + log10
mSet = Normalization(mSet, "SumNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

# sanity check
mSet = SanityCheckData(mSet)

# univariate analysis
mSet = PlotNormSummary(mSet, "univariate_norm_3_", "png", 72, width=NA)
mSet = PlotSampleNormSummary(mSet, "univariate_snorm_3_", "png", 72, width=NA)
# t-test + FDR
# mSet = Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
# mSet = PlotTT(mSet, "univariate_tt_0_", "png", 72, width=NA)
# non-parametric test + FDR
mSet = Ttests.Anal(mSet, T, 0.10, FALSE, TRUE, "fdr", FALSE)
mSet = PlotTT(mSet, "univariate_tt_1_", "png", 72, width=NA)
# non-parametric test + RAW (without FDR)
# mSet = Ttests.Anal(mSet, T, 0.05, FALSE, TRUE, "raw", FALSE)
# mSet = PlotTT(mSet, "univariate_tt_2_", "png", 72, width=NA)
# non-parametric test + RAW (without FDR) + unequal variance
# mSet = Ttests.Anal(mSet, T, 0.05, FALSE, FALSE, "raw", FALSE)
# mSet = PlotTT(mSet, "univariate_tt_3_", "png", 72, width=NA)
