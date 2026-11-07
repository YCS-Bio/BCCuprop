#!/usr/bin/env Rscript
# Complete Immune Reprogramming Analysis for 5 Genes
# No external source files needed

# ============================================================================
# 1. LOAD LIBRARIES
# ============================================================================

library(tidyverse)
library(ggplot2)
library(reshape2)
library(corrplot)
library(survival)
library(survminer)
library(ggpubr)

# Load deconvolution libraries
library(immunedeconv)  # Main package for multiple methods
library(xCell)         # For xCell specifically

# ============================================================================
# 2. LOAD AND PREPARE YOUR DATA
# ============================================================================

cat("Step 1: Loading your expression data...\n")

# CHANGE THIS to your actual file path
expression_file <- "your_expression_data.txt"

# Read your data (adjust parameters based on your file format)
# Assuming: tab-delimited, genes as rows, samples as columns
# First column contains gene symbols
expression_data <- read.delim(expression_file, 
                              header = TRUE, 
                              row.names = 1,  # First column as row names
                              check.names = FALSE)  # Preserve column names as-is

# Convert to matrix
exp_matrix <- as.matrix(expression_data)

# Check dimensions
cat(sprintf("Expression matrix dimensions: %d genes x %d samples\n", 
            nrow(exp_matrix), ncol(exp_matrix)))

# ============================================================================
# 3. DEFINE YOUR 5 GENES
# ============================================================================

# REPLACE THESE WITH YOUR ACTUAL GENE SYMBOLS
my_genes <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")

# Check which genes are in your data
available_genes <- my_genes[my_genes %in% rownames(exp_matrix)]
cat(sprintf("Found %d of %d genes in your data: %s\n", 
            length(available_genes), 
            length(my_genes), 
            paste(available_genes, collapse = ", ")))

if (length(available_genes) == 0) {
  stop("None of your genes were found in the expression data. Please check gene symbols.")
}

# Extract expression for your genes
gene_exp <- exp_matrix[available_genes, , drop = FALSE]

# ============================================================================
# 4. IMMUNE CELL DECONVOLUTION (using immunedeconv package)
# ============================================================================

cat("\nStep 2: Performing immune cell deconvolution...\n")

# Method 1: CIBERSORT (requires internet connection for first run)
cat("Running CIBERSORT...\n")
tryCatch({
  cibersort_scores <- deconvolute(exp_matrix, "cibersort", arrays = FALSE)
}, error = function(e) {
  cat("CIBERSORT failed:", e$message, "\n")
  cibersort_scores <- NULL
})

# Method 2: xCell
cat("Running xCell...\n")
tryCatch({
  xcell_scores <- deconvolute(exp_matrix, "xcell")
}, error = function(e) {
  cat("xCell failed:", e$message, "\n")
  # Alternative: use xCell package directly
  xcell_scores <- as.data.frame(xCellAnalysis(exp_matrix))
  xcell_scores <- as.data.frame(t(xcell_scores))
})

# Method 3: EPIC
cat("Running EPIC...\n")
tryCatch({
  epic_scores <- deconvolute(exp_matrix, "epic")
}, error = function(e) {
  cat("EPIC failed:", e$message, "\n")
  epic_scores <- NULL
})

# Method 4: MCPcounter (built into immunedeconv)
cat("Running MCPcounter...\n")
tryCatch({
  mcpcounter_scores <- deconvolute(exp_matrix, "mcp_counter")
}, error = function(e) {
  cat("MCPcounter failed:", e$message, "\n")
  mcpcounter_scores <- NULL
})

# ============================================================================
# 5. COMBINE AND PROCESS IMMUNE SCORES
# ============================================================================

cat("\nStep 3: Combining immune scores...\n")

# Function to process deconvolution results
process_scores <- function(scores, prefix = "") {
  if (is.null(scores)) return(NULL)
  
  # Transpose if needed (different methods have different formats)
  if ("cell_type" %in% colnames(scores)) {
    # immunedeconv format: cell_type in first column
    cell_types <- scores$cell_type
    scores_mat <- as.matrix(scores[, -1, drop = FALSE])
    rownames(scores_mat) <- cell_types
    scores_df <- as.data.frame(t(scores_mat))
  } else {
    # Already samples as rows, cell types as columns
    scores_df <- as.data.frame(scores)
  }
  
  # Add prefix to column names
  if (prefix != "") {
    colnames(scores_df) <- paste0(prefix, "_", colnames(scores_df))
  }
  
  return(scores_df)
}

# Process each method
xcell_processed <- process_scores(xcell_scores, "xCell")
cibersort_processed <- process_scores(cibersort_scores, "CIBERSORT")
epic_processed <- process_scores(epic_scores, "EPIC")
mcpcounter_processed <- process_scores(mcpcounter_scores, "MCPcounter")

# Combine all scores
immune_scores <- data.frame()

# Start with xCell (most comprehensive)
if (!is.null(xcell_processed)) {
  immune_scores <- xcell_processed
}

# Add other methods
if (!is.null(cibersort_processed)) {
  common_samples <- intersect(rownames(immune_scores), rownames(cibersort_processed))
  if (length(common_samples) > 0) {
    immune_scores <- cbind(immune_scores[common_samples, , drop = FALSE],
                           cibersort_processed[common_samples, , drop = FALSE])
  }
}

# Select key immune cell types for breast cancer analysis
key_cell_types <- c(
  # T cells
  "xCell_CD8+ T-cells", "xCell_CD4+ T-cells", "xCell_Tregs",
  "CIBERSORT_T cells CD8", "CIBERSORT_T cells regulatory (Tregs)",
  
  # Macrophages
  "xCell_Macrophages", "xCell_Macrophages M1", "xCell_Macrophages M2",
  "CIBERSORT_Macrophages M0", "CIBERSORT_Macrophages M1", "CIBERSORT_Macrophages M2",
  
  # Other immune cells
  "xCell_B-cells", "xCell_NK cells", "xCell_Neutrophils",
  "xCell_Dendritic cells", "xCell_Mast cells",
  
  # Stromal cells
  "xCell_Cancer associated fibroblasts", "xCell_Endothelial cells",
  
  # MCPcounter if available
  "MCPcounter_CD8 T cells", "MCPcounter_Cytotoxic lymphocytes",
  "MCPcounter_Monocytic lineage", "MCPcounter_Fibroblasts"
)

# Filter to available cell types
available_cells <- key_cell_types[key_cell_types %in% colnames(immune_scores)]
immune_scores <- immune_scores[, available_cells, drop = FALSE]

# Save immune scores
write.csv(immune_scores, "immune_cell_scores.csv", row.names = TRUE)
cat(sprintf("Saved immune scores for %d cell types across %d samples\n", 
            ncol(immune_scores), nrow(immune_scores)))

# ============================================================================
# 6. CORRELATION ANALYSIS
# ============================================================================

cat("\nStep 4: Calculating gene-immune cell correlations...\n")

# Prepare data for correlation analysis
# Your gene expression (genes as rows, samples as columns)
gene_exp_t <- as.data.frame(t(gene_exp))  # Transpose: samples as rows

# Ensure same samples in both datasets
common_samples <- intersect(rownames(gene_exp_t), rownames(immune_scores))
cat(sprintf("Using %d common samples for analysis\n", length(common_samples)))

if (length(common_samples) < 10) {
  warning("Very few common samples! Check sample names in both datasets.")
}

gene_exp_t <- gene_exp_t[common_samples, , drop = FALSE]
immune_scores <- immune_scores[common_samples, , drop = FALSE]

# Calculate Spearman correlations
cor_results <- data.frame()

for (gene in colnames(gene_exp_t)) {
  for (cell_type in colnames(immune_scores)) {
    # Calculate correlation
    cor_test <- cor.test(gene_exp_t[[gene]], 
                         immune_scores[[cell_type]], 
                         method = "spearman",
                         exact = FALSE)
    
    # Store results
    cor_results <- rbind(cor_results, data.frame(
      Gene = gene,
      Immune_Cell = cell_type,
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# Adjust p-values for multiple testing
cor_results$FDR <- p.adjust(cor_results$P_value, method = "BH")

# Sort by absolute correlation
cor_results <- cor_results[order(-abs(cor_results$Correlation)), ]

# Save results
write.csv(cor_results, "gene_immune_correlations.csv", row.names = FALSE)
cat(sprintf("Calculated %d correlations. Saved to file.\n", nrow(cor_results)))

# ============================================================================
# 7. VISUALIZATION
# ============================================================================

cat("\nStep 5: Creating visualizations...\n")

# 7.1 Heatmap of top correlations
pdf("correlation_heatmap.pdf", width = 12, height = 10)

# Select top 20 correlations
top_n <- min(20, nrow(cor_results))
top_cors <- cor_results[1:top_n, ]

# Create matrix for heatmap
heatmap_data <- dcast(top_cors, Gene ~ Immune_Cell, value.var = "Correlation")
rownames(heatmap_data) <- heatmap_data$Gene
heatmap_matrix <- as.matrix(heatmap_data[, -1])

# Create heatmap
corrplot(heatmap_matrix,
         method = "color",
         type = "full",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 45,
         mar = c(0, 0, 2, 0),
         title = paste("Top", top_n, "Gene-Immune Cell Correlations"),
         addCoef.col = "black",
         number.cex = 0.7,
         cl.cex = 0.8,
         is.corr = FALSE)

dev.off()

# 7.2 Scatter plots for top 4 correlations
pdf("top_correlations_scatter.pdf", width = 12, height = 10)

top_4 <- head(cor_results, 4)
par(mfrow = c(2, 2))

for (i in 1:nrow(top_4)) {
  gene <- top_4$Gene[i]
  cell <- top_4$Immune_Cell[i]
  cor_val <- round(top_4$Correlation[i], 3)
  p_val <- formatC(top_4$P_value[i], format = "e", digits = 2)
  
  plot(gene_exp_t[[gene]], immune_scores[[cell]],
       xlab = paste(gene, "Expression"),
       ylab = paste(cell, "Score"),
       main = paste(gene, "vs", cell),
       sub = paste("r =", cor_val, "p =", p_val),
       pch = 19, col = rgb(0.2, 0.4, 0.8, 0.6),
       cex = 1.2)
  
  # Add regression line
  abline(lm(immune_scores[[cell]] ~ gene_exp_t[[gene]]), 
         col = "red", lwd = 2)
}

dev.off()

# 7.3 Bar plot of strongest correlations per gene
pdf("gene_immune_summary.pdf", width = 10, height = 6)

# Get strongest positive and negative correlation for each gene
gene_summary <- cor_results %>%
  group_by(Gene) %>%
  summarize(
    Max_Positive = max(Correlation, na.rm = TRUE),
    Max_Negative = min(Correlation, na.rm = TRUE),
    Pos_Cell = Immune_Cell[which.max(Correlation)],
    Neg_Cell = Immune_Cell[which.min(Correlation)]
  ) %>%
  pivot_longer(cols = c(Max_Positive, Max_Negative),
               names_to = "Direction",
               values_to = "Correlation") %>%
  mutate(Cell_Type = ifelse(Direction == "Max_Positive", Pos_Cell, Neg_Cell))

# Create bar plot
ggplot(gene_summary, aes(x = Gene, y = Correlation, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = Cell_Type, y = Correlation + sign(Correlation) * 0.02),
            position = position_dodge(width = 0.7),
            size = 3, angle = 45, hjust = 0) +
  labs(title = "Strongest Immune Correlations for Each Gene",
       y = "Correlation Coefficient",
       x = "") +
  scale_fill_manual(values = c("Max_Positive" = "red3", "Max_Negative" = "blue3"),
                    labels = c("Most Positive", "Most Negative")) +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold", size = 11, angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"))

dev.off()

# 7.4 Network-style visualization
pdf("correlation_network.pdf", width = 14, height = 10)

# Filter significant correlations (FDR < 0.1)
sig_cors <- cor_results[cor_results$FDR < 0.1, ]

if (nrow(sig_cors) > 0) {
  # Create node data
  genes <- unique(sig_cors$Gene)
  cells <- unique(sig_cors$Immune_Cell)
  
  # Create edge data
  edges <- sig_cors %>%
    mutate(width = abs(Correlation) * 3,
           color = ifelse(Correlation > 0, "red", "blue"),
           linetype = ifelse(Correlation > 0, "solid", "dashed"))
  
  # Simple visualization using segments
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, length(genes) + length(cells) + 2),
       xlab = "", ylab = "", axes = FALSE, main = "Significant Gene-Immune Cell Correlations")
  
  # Place genes on left, immune cells on right
  gene_y <- seq(2, by = 2, length.out = length(genes))
  cell_y <- seq(2, by = 2, length.out = length(cells))
  
  # Draw genes (left)
  points(rep(2, length(genes)), gene_y, pch = 19, col = "darkgreen", cex = 3)
  text(rep(1.8, length(genes)), gene_y, genes, adj = 1, cex = 0.8, font = 2)
  
  # Draw immune cells (right)
  points(rep(8, length(cells)), cell_y, pch = 15, col = "purple", cex = 3)
  text(rep(8.2, length(cells)), cell_y, cells, adj = 0, cex = 0.7)
  
  # Draw edges for significant correlations
  for (i in 1:nrow(edges)) {
    gene_idx <- which(genes == edges$Gene[i])
    cell_idx <- which(cells == edges$Immune_Cell[i])
    
    lines(c(2, 8), c(gene_y[gene_idx], cell_y[cell_idx]),
          col = edges$color[i],
          lwd = edges$width[i],
          lty = ifelse(edges$Correlation[i] > 0, "solid", "dashed"))
  }
  
  # Add legend
  legend("topright", 
         legend = c("Positive correlation", "Negative correlation", "Gene", "Immune Cell"),
         col = c("red", "blue", "darkgreen", "purple"),
         pch = c(NA, NA, 19, 15),
         lty = c("solid", "dashed", NA, NA),
         lwd = 2)
} else {
  plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No significant correlations (FDR < 0.1)", cex = 1.2)
}

dev.off()

# ============================================================================
# 8. CREATE SUMMARY STATISTICS AND REPORT
# ============================================================================

cat("\nStep 6: Generating summary report...\n")

# Calculate summary statistics
summary_stats <- cor_results %>%
  group_by(Gene) %>%
  summarize(
    Num_Significant_Pos = sum(Correlation > 0.3 & FDR < 0.1, na.rm = TRUE),
    Num_Significant_Neg = sum(Correlation < -0.3 & FDR < 0.1, na.rm = TRUE),
    Mean_Correlation = mean(Correlation, na.rm = TRUE),
    Strongest_Positive = max(Correlation, na.rm = TRUE),
    Strongest_Negative = min(Correlation, na.rm = TRUE),
    Most_Positive_Cell = Immune_Cell[which.max(Correlation)],
    Most_Negative_Cell = Immune_Cell[which.min(Correlation)]
  )

# Save summary
write.csv(summary_stats, "correlation_summary_statistics.csv", row.names = FALSE)

# Generate text report
sink("analysis_report.txt")

cat("=", strrep("=", 70), "\n", sep = "")
cat("IMMUNE REPROGRAMMING ANALYSIS REPORT\n")
cat("=", strrep("=", 70), "\n\n")

cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("Time:", format(Sys.time(), "%H:%M:%S"), "\n\n")

cat("ANALYSIS PARAMETERS\n")
cat(strrep("-", 70), "\n")
cat("Genes analyzed:", paste(available_genes, collapse = ", "), "\n")
cat("Samples analyzed:", nrow(immune_scores), "\n")
cat("Immune cell types:", ncol(immune_scores), "\n")
cat("Deconvolution methods used:", 
    paste(c("xCell", "CIBERSORT", "EPIC", "MCPcounter")[
      c(!is.null(xcell_scores), !is.null(cibersort_scores), 
        !is.null(epic_scores), !is.null(mcpcounter_scores))
    ], collapse = ", "), "\n\n")

cat("SUMMARY STATISTICS\n")
cat(strrep("-", 70), "\n")
print(summary_stats)

cat("\nTOP SIGNIFICANT CORRELATIONS (FDR < 0.1)\n")
cat(strrep("-", 70), "\n")
sig_top <- cor_results[cor_results$FDR < 0.1, ]
if (nrow(sig_top) > 0) {
  print(head(sig_top[order(-abs(sig_top$Correlation)), ], 10))
} else {
  cat("No significant correlations at FDR < 0.1\n")
}

cat("\nIMMUNOSUPPRESSIVE ASSOCIATIONS\n")
cat(strrep("-", 70), "\n")

# Check associations with immunosuppressive cells
immunosuppressive_cells <- c("Macrophages M2", "Tregs", "T cells regulatory", 
                             "Cancer associated fibroblasts", "Myeloid derived suppressor")

for (gene in available_genes) {
  gene_cors <- cor_results[cor_results$Gene == gene & cor_results$FDR < 0.1, ]
  
  # Find immunosuppressive associations
  immunosuppressive <- gene_cors[grepl(paste(immunosuppressive_cells, collapse = "|"), 
                                       gene_cors$Immune_Cell, ignore.case = TRUE), ]
  
  if (nrow(immunosuppressive) > 0) {
    cat(sprintf("\n%s shows immunosuppressive associations:\n", gene))
    for (i in 1:nrow(immunosuppressive)) {
      cat(sprintf("  • %s: r = %.2f, FDR = %.2e\n", 
                  immunosuppressive$Immune_Cell[i],
                  immunosuppressive$Correlation[i],
                  immunosuppressive$FDR[i]))
    }
  }
}

cat("\n" + strrep("=", 70) + "\n")
cat("INTERPRETATION GUIDELINES:\n")
cat("• Positive correlation with M2 macrophages/Tregs = Possible immunosuppressive role\n")
cat("• Negative correlation with CD8+ T cells/NK cells = Possible immune evasion role\n")
cat("• Strong correlations (>|0.4|) are biologically meaningful\n")
cat("• Check individual scatter plots for non-linear relationships\n")
cat(strrep("=", 70) + "\n")

sink()

# Also print summary to console
cat("\nANALYSIS COMPLETE!\n")
cat("Generated files:\n")
cat("1. immune_cell_scores.csv - Immune cell abundance for each sample\n")
cat("2. gene_immune_correlations.csv - All correlation results\n")
cat("3. correlation_heatmap.pdf - Heatmap of top correlations\n")
cat("4. top_correlations_scatter.pdf - Scatter plots for top correlations\n")
cat("5. gene_immune_summary.pdf - Bar plot summary\n")
cat("6. correlation_network.pdf - Network visualization\n")
cat("7. correlation_summary_statistics.csv - Summary statistics\n")
cat("8. analysis_report.txt - Detailed analysis report\n\n")

# Print key findings
cat("KEY FINDINGS:\n")
for (i in 1:min(5, nrow(cor_results))) {
  row <- cor_results[i, ]
  cat(sprintf("%s vs %s: r = %.2f (FDR = %.2e)\n", 
              row$Gene, row$Immune_Cell, row$Correlation, row$FDR))
}