expr=read.table("TCGA_TPM.txt",row.names = 1, check.names = F, sep = "\t",header=T )
y=read.table("y.txt",row.names = 1, check.names = F, sep = "\t",header=T )
expr=expr[rownames(y)]
y <- y[colnames(expr), 1]   # 按 expr 列名排序，取第一列
y <- factor(y)                  # 转为因子（适合分类分析）

# expr = TPM matrix (genes x samples)

# log-transform TPM
log_expr <- log2(expr + 1)

# compute variance per gene
gene_var <- apply(log_expr, 1, var)

# keep only genes with variance > 0
log_expr_filt <- log_expr[gene_var > 0, ]





library(ggplot2)

run_pca <- function(expr, y) {
  # log-transform TPM to stabilize variance
  log_expr <- log2(expr + 1)
  
  # transpose to samples × genes
  pca <- prcomp(t(log_expr), scale. = TRUE)
  
  df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    y = factor(y)
  )
  
  ggplot(df, aes(PC1, PC2, color = y)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)"),
      color = "Class"
    )
}
library(umap)
library(ggplot2)

run_umap <- function(expr, y, n_neighbors = 15, min_dist = 0.1) {
  log_expr <- log2(expr + 1)
  
  um <- umap(t(log_expr), n_neighbors = n_neighbors, min_dist = min_dist)
  
  df <- data.frame(
    UMAP1 = um$layout[, 1],
    UMAP2 = um$layout[, 2],
    y = factor(y)
  )
  
  ggplot(df, aes(UMAP1, UMAP2, color = y)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(color = "Class")
}
p1 <- run_pca(log_expr_filt, y)
p2 <- run_umap(log_expr_filt, y)

p1
p2

pca <- prcomp(t(log_expr_filt), scale. = TRUE)

df_pca <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  y = factor(y)
)

ggplot(df_pca, aes(PC1, PC2, color = y)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14)

head(pca$rotation[order(abs(pca$rotation[,1]), decreasing=TRUE), 1], 30)
head(pca$rotation[order(abs(pca$rotation[,2]), decreasing=TRUE), 2], 30)
