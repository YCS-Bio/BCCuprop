# Install if needed
if (!requireNamespace("ConsensusClusterPlus", quietly=TRUE)) {
  BiocManager::install("ConsensusClusterPlus")
}
library(ConsensusClusterPlus)

# 1 Load data (columns = samples)
expr <- read.delim("TCGA_TPM_patient_1097_train.txt", check.names=FALSE)
gene_col <- colnames(expr)[1]
mat <- as.matrix(expr[ , -1])        # remove gene id column
rownames(mat) <- expr[[1]]

# 2 Preprocess
mat <- log2(mat + 1)                 # transform
# optional: filter top variable genes
vars <- apply(mat, 1, var)
mat <- mat[order(vars, decreasing=TRUE)[1:5000], ]

# 3 Run ConsensusClusterPlus
results <- ConsensusClusterPlus(
  d = mat,
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  innerLinkage = "average",
  finalLinkage = "average",
  distance = "pearson",
  seed = 1234,
  title = "consensus_patient_1097",
  plot = "png",
  writeTable = TRUE
)