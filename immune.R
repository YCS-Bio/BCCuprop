# ---------------------------
# ssGSEA (C7/C8) + xCell pipeline
# ---------------------------

# Environment
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())
setwd("D:\\1-7cuproptosisenrichment")

# Install / load packages (run once)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("data.table","GSVA","GSEABase","msigdbr","limma","xCell","pheatmap","ggpubr","reshape2","dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  library(p, character.only = TRUE)
}

# ---------------------------
# 1. Load expression and cluster files
# ---------------------------
expr_file <- "TCGA_TPM.txt"                # edit if needed
cluster_file <- "Cuproptosiscluster2.txt"

# Read expression: first column gene symbol, rest samples
rt <- fread(expr_file, data.table = FALSE)
rownames(rt) <- rt[,1]
expr <- as.matrix(rt[,-1])
mode(expr) <- "numeric"

# Average duplicate gene symbols
expr <- avereps(expr)

# Optional: if raw counts, transform
expr <- log2(expr + 0.1)

# Read cluster labels (row names must be sample IDs)
cluster <- read.table(cluster_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# Ensure cluster column name
if (!"Cuproptosiscluster" %in% colnames(cluster)) {
  colnames(cluster)[1] <- "Cuproptosiscluster"
}

# ---------------------------
# 2. Download C7 and C8 via msigdbr and build gene lists
# ---------------------------
message("Downloading MSigDB C7 and C8 via msigdbr...")
c7_df <- msigdbr(species = "Homo sapiens", category = "C7")
c8_df <- msigdbr(species = "Homo sapiens", category = "C8")

# Convert to named list: gs_name -> vector of gene symbols
c7_list <- split(c7_df$gene_symbol, c7_df$gs_name)
c8_list <- split(c8_df$gene_symbol, c8_df$gs_name)

# ---------------------------
# 3. Prune gene sets and align to expression genes
# ---------------------------
prune_gs <- function(glist, genes, min_sz = 15, max_sz = 500) {
  glist2 <- lapply(glist, function(x) unique(intersect(x, genes)))
  sizes <- sapply(glist2, length)
  keep <- which(sizes >= min_sz & sizes <= max_sz)
  glist2[keep]
}

genes_available <- rownames(expr)
c7_pruned <- prune_gs(c7_list, genes_available, min_sz = 15, max_sz = 500)
c8_pruned <- prune_gs(c8_list, genes_available, min_sz = 10, max_sz = 500)  # C8 sets can be smaller

message("C7 sets kept: ", length(c7_pruned))
message("C8 sets kept: ", length(c8_pruned))

# ---------------------------
# 4. Run ssGSEA (GSVA implementation)
# ---------------------------
set.seed(123)
# Use kcdf = "Gaussian" for log-like continuous data; use "Poisson" for raw counts
scores_c7 <- gsva(expr, c7_pruned, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, ssgsea.norm = TRUE)
scores_c8 <- gsva(expr, c8_pruned, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, ssgsea.norm = TRUE)

# Save raw ssGSEA scores
write.table(rbind(id = colnames(scores_c7), scores_c7), file = "ssGSEA_C7_scores.txt", sep = "\t", quote = FALSE, col.names = FALSE)
write.table(rbind(id = colnames(scores_c8), scores_c8), file = "ssGSEA_C8_scores.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# ---------------------------
# 5. Run xCell for complementary cell-type enrichment
# ---------------------------
# xCell expects genes as rows and samples as columns
message("Running xCell (this may take a minute)...")
xcell_res <- xCellAnalysis(expr)   # returns matrix: cell types x samples
write.table(rbind(id = colnames(xcell_res), xcell_res), file = "xCell_scores.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# ---------------------------
# 6. Prepare merged data for comparison with clusters
# ---------------------------
# Transpose ssGSEA results to samples x geneSets
ssgsea_c7_t <- t(scores_c7)
ssgsea_c8_t <- t(scores_c8)
xcell_t <- t(xcell_res)

# Intersect samples with cluster labels
common_samples <- Reduce(intersect, list(rownames(ssgsea_c7_t), rownames(ssgsea_c8_t), rownames(xcell_t), rownames(cluster)))
ssgsea_c7_t <- ssgsea_c7_t[common_samples, , drop = FALSE]
ssgsea_c8_t <- ssgsea_c8_t[common_samples, , drop = FALSE]
xcell_t <- xcell_t[common_samples, , drop = FALSE]
cluster_sub <- cluster[common_samples, , drop = FALSE]

# Combine into data frames for plotting
ssgsea_c7_df <- as.data.frame(ssgsea_c7_t); ssgsea_c7_df$Sample <- rownames(ssgsea_c7_df)
ssgsea_c7_df <- left_join(ssgsea_c7_df, data.frame(Sample = rownames(cluster_sub), cluster_sub), by = "Sample")

ssgsea_c8_df <- as.data.frame(ssgsea_c8_t); ssgsea_c8_df$Sample <- rownames(ssgsea_c8_df)
ssgsea_c8_df <- left_join(ssgsea_c8_df, data.frame(Sample = rownames(cluster_sub), cluster_sub), by = "Sample")

xcell_df <- as.data.frame(xcell_t); xcell_df$Sample <- rownames(xcell_df)
xcell_df <- left_join(xcell_df, data.frame(Sample = rownames(cluster_sub), cluster_sub), by = "Sample")

# ---------------------------
# 7. Statistical comparison between two clusters
# ---------------------------
# Function to compute Wilcoxon p-values and adjusted p-values for each feature
compare_groups <- function(mat, group_vec, test = "wilcox") {
  features <- colnames(mat)
  pvals <- sapply(features, function(f) {
    x <- mat[, f]
    g1 <- x[group_vec == unique(group_vec)[1]]
    g2 <- x[group_vec == unique(group_vec)[2]]
    if (length(unique(group_vec)) != 2) return(NA)
    res <- wilcox.test(g1, g2)$p.value
    res
  })
  p.adj <- p.adjust(pvals, method = "BH")
  data.frame(feature = features, p.value = pvals, p.adj = p.adj, stringsAsFactors = FALSE)
}

group_vec <- cluster_sub$Cuproptosiscluster

# For speed, test a subset or all features
res_c7 <- compare_groups(ssgsea_c7_t, group_vec)
res_c8 <- compare_groups(ssgsea_c8_t, group_vec)
res_xcell <- compare_groups(xcell_t, group_vec)

# Save results
write.table(res_c7, file = "ssGSEA_C7_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res_c8, file = "ssGSEA_C8_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res_xcell, file = "xCell_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 8. Visualization examples
# ---------------------------

# 8A. Heatmap of top variable gene sets (C7)
top_n <- 30
var_order <- order(apply(ssgsea_c7_t, 2, var), decreasing = TRUE)
top_sets <- colnames(ssgsea_c7_t)[var_order[1:top_n]]
ann_col <- data.frame(Cluster = cluster_sub$Cuproptosiscluster)
rownames(ann_col) <- rownames(cluster_sub)

pdf("heatmap_ssGSEA_C7_top30.pdf", width = 10, height = 8)
pheatmap(t(ssgsea_c7_t[, top_sets]), scale = "row", annotation_col = ann_col, show_colnames = FALSE, fontsize_row = 7)
dev.off()

# 8B. Boxplots for selected features side-by-side (C8 and xCell)
# Choose a few biologically relevant sets/cell types to plot
select_c8 <- head(order(apply(ssgsea_c8_t, 2, function(x) abs(mean(x[group_vec==unique(group_vec)[1]]) - mean(x[group_vec==unique(group_vec)[2]]))), decreasing = TRUE), 8)
select_xcell <- head(order(apply(xcell_t, 2, function(x) abs(mean(x[group_vec==unique(group_vec)[1]]) - mean(x[group_vec==unique(group_vec)[2]]))), decreasing = TRUE), 8)

plot_boxsets <- function(df, features, cluster_col = "Cuproptosiscluster", out_pdf = "boxplot.pdf", title = "") {
  long <- melt(df[, c(features, cluster_col)], id.vars = cluster_col)
  colnames(long) <- c("Cluster", "Feature", "Score")
  long$Feature <- factor(long$Feature, levels = unique(long$Feature))
  p <- ggboxplot(long, x = "Feature", y = "Score", color = "Cluster", add = "jitter", palette = "jco") +
    rotate_x_text(45) +
    stat_compare_means(aes(group = Cluster), method = "wilcox.test", label = "p.signif")
  pdf(out_pdf, width = 12, height = 6)
  print(p + ggtitle(title))
  dev.off()
}

plot_boxsets(ssgsea_c8_df, colnames(ssgsea_c8_df)[select_c8], out_pdf = "boxplot_ssGSEA_C8_selected.pdf", title = "Top differential C8 sets")
plot_boxsets(xcell_df, colnames(xcell_df)[select_xcell], out_pdf = "boxplot_xCell_selected.pdf", title = "Top differential xCell cell types")

# ---------------------------
# 9. Quick concordance check between ssGSEA cell-type sets and xCell
# ---------------------------
# Example: correlate a C8 cell-type set with xCell estimate for same cell type
# Build correlation table for matching names (best-effort string matching)
cor_table <- data.frame()
for (c8_name in colnames(ssgsea_c8_t)) {
  for (xc_name in colnames(xcell_t)) {
    cor_val <- cor(ssgsea_c8_t[, c8_name], xcell_t[, xc_name], method = "spearman")
    cor_table <- rbind(cor_table, data.frame(C8 = c8_name, xCell = xc_name, spearman = cor_val))
  }
}
# Keep top correlations
top_cor <- cor_table[order(abs(cor_table$spearman), decreasing = TRUE), ][1:50, ]
write.table(top_cor, file = "C8_xCell_correlations_top50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# End of script
# ---------------------------
message("Pipeline finished. Outputs: ssGSEA_C7_scores.txt, ssGSEA_C8_scores.txt, xCell_scores.txt, stats files, heatmaps and boxplots.")