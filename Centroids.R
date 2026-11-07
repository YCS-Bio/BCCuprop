#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
setwd("D:/ML_REO/split0221")
workDir <- 'D:/ML_REO/split0221'

#=======================================================
# 1.Build Centroids From Training Data
#=======================================================
# Load TCGA 19-gene expression (genes x samples)
expr_tcga <- read.table("expr_crg19_train_scaled.txt", sep="\t", header=T, row.names=1, check.names=FALSE)

# Load cluster labels from ConsensusClusterPlus
clusters <- read.table("CuproptosisCluster4.txt", sep="\t", header=TRUE)  # sample_id, cluster

# Reorder expression to match cluster file
expr_tcga <- expr_tcga[, clusters$ID]

# Compute centroids (mean expression per cluster)

# Determine number of clusters automatically
K <- length(unique(clusters$Cuproptosiscluster))

# Compute centroids
centroids <- sapply(1:K, function(k) {
  rowMeans(expr_tcga[, clusters$Cuproptosiscluster == k, drop=FALSE])
})

colnames(centroids) <- paste0("Cluster", 1:K)


write.table(centroids, "CRG19_k4_centroids.txt", sep="\t", quote=FALSE)

#========================================================
# 2. Prepare New Dataset (RNA‑seq or Microarray)
#========================================================
expr <- read.delim("expr_crg19_test.txt", check.names=FALSE)
gene_col <- colnames(expr)[1]
mat <- as.matrix(expr[ , -1])        # remove gene id column
rownames(mat) <- expr[[1]]

# 2 Preprocess
mat <- log2(mat + 1)                 # transform
# optional: filter top variable genes
vars <- apply(mat, 1, var)
mat <- mat[order(vars, decreasing=TRUE)[1:19], ]

#  z score by gene using scaler from training set

# 加载参数
scaler <- readRDS("scaler.rds")
center <- scaler$center
scale  <- scaler$scale

# 新数据
new <- t(mat)

# 获取新数据的列名
new_cols <- colnames(new)

# 确保所有列在训练集中存在（否则报错或处理缺失）
if (!all(new_cols %in% names(center))) {
  stop("新数据包含训练集中未出现的列")
}

# 按新数据列名重新排列 center 和 scale
center_aligned <- center[new_cols]
scale_aligned  <- scale[new_cols]

# 应用标准化（现在顺序一致）
scaled_new <- scale(new, center = center_aligned, scale = scale_aligned)

mat_z <- t(scaled_new)
write.table(mat_z, file="expr_crg19_test_scaled.txt", sep="\t", quote=F, col.names=T)



expr_new <- read.table("expr_crg19_test_scaled.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)

# Ensure gene order matches centroid
centroids=read.table("CRG19_k4_centroids.txt", sep="\t", row.names=1)
expr_new <- expr_new[rownames(centroids), , drop=FALSE]
#=====================================================
# 3. Assign New Samples to Clusters
#=====================================================
assign_cluster <- function(sample_vec, centroids, method="pearson") {
  sims <- apply(centroids, 2, function(cent) cor(sample_vec, cent, method=method))
  best <- names(which.max(sims))
  return(list(cluster=best, similarity=sims))
}
assignments <- sapply(colnames(expr_new), function(s) {
  assign_cluster(expr_new[,s], centroids, method="pearson")$cluster
})

write.table(assignments, "NEW_DATASET_cluster_assignments3.txt", sep="\t", quote=FALSE)
