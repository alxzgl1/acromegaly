# LINK: https://www.metaboanalyst.ca/MetaboAnalyst/docs/RTutorial.xhtml

# load package
library("MetaboAnalystR")

# set directory
setwd("D:/code/__R/Seedcorn") # L15 laptop

# filename
aFilename = "D:/data/acromegaly/import/Lipids_P_S30F25_NA.csv" # HILIC_N, HILIC_P, Lipids_N, Lipids_P

# init
mSet = InitDataObjects("pktable", "stat", FALSE)
mSet = Read.TextData(mSet, aFilename, "colu", "disc") 

# univariate / multivariate
bUnivariate = FALSE
if (bUnivariate == TRUE) {
  print("Univariate analysis", quote=FALSE)
  # imputation with min
  mSet = ImputeMissingVar(mSet, method="min")
  # prenorm
  mSet = PreparePrenormData(mSet)
  # normalise
  mSet = Normalization(mSet, "SumNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
  # save
  mSet = SanityCheckData(mSet)
  SaveTransformedData(mSet) 
} else { 
  print("Multivariate analysis", quote=FALSE)
  # imputation with KNN
  mSet = ImputeMissingVar(mSet, method="knn_var")
  # prenorm
  mSet = PreparePrenormData(mSet)
  # normalise
  mSet = Normalization(mSet, "SumNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
  # save
  mSet = SanityCheckData(mSet)
  SaveTransformedData(mSet) 
}

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
