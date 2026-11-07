#!/usr/bin/env Rscript
# Complete Immune Reprogramming Analysis for 5 Genes
# Only uses freely available methods (no CIBERSORT required)

# ============================================================================
# 1. LOAD LIBRARIES (ONLY FREE METHODS)
# ============================================================================

library(tidyverse)
library(ggplot2)
library(reshape2)
library(corrplot)
library(survival)
library(survminer)
library(ggpubr)

# Load deconvolution libraries for free methods
library(immunedeconv)  # For EPIC, MCPcounter, quantiseq
library(xCell)         # For xCell

# ============================================================================
# 2. LOAD YOUR DATA
# ============================================================================

cat("Step 1: Loading your expression data...\n")

# CHANGE THIS to your actual file path
expression_file <- "TCGA_TPM.txt"

# Read your data
expression_data <- read.delim(expression_file, 
                              header = TRUE, 
                              row.names = 1,
                              check.names = FALSE)

# Convert to matrix
exp_matrix <- as.matrix(expression_data)

# Check dimensions
cat(sprintf("Expression matrix dimensions: %d genes x %d samples\n", 
            nrow(exp_matrix), ncol(exp_matrix)))

# Check data range
cat(sprintf("Expression range: %.2f to %.2f\n", 
            min(exp_matrix, na.rm = TRUE), 
            max(exp_matrix, na.rm = TRUE)))

# ============================================================================
# 3. DEFINE YOUR 5 GENES
# ============================================================================

# REPLACE THESE WITH YOUR ACTUAL GENE SYMBOLS
my_genes <- c("NIPBL", "SEC24C", "PRDM10", "LTN1", "IREB2")

# Check which genes are in your data
available_genes <- my_genes[my_genes %in% rownames(exp_matrix)]
cat(sprintf("Found %d of %d genes in your data: %s\n", 
            length(available_genes), 
            length(my_genes), 
            paste(available_genes, collapse = ", ")))

if (length(available_genes) == 0) {
  cat("\nAvailable gene symbols (first 20):\n")
  print(head(rownames(exp_matrix), 20))
  stop("None of your genes were found in the expression data.")
}

# Extract expression for your genes
gene_exp <- exp_matrix[available_genes, , drop = FALSE]

# ============================================================================
# 4. IMMUNE CELL DECONVOLUTION (FREE METHODS ONLY)
# ============================================================================

cat("\nStep 2: Performing immune cell deconvolution...\n")

# Function to run deconvolution with error handling
safe_deconvolute <- function(method_name, exp_matrix) {
  cat(sprintf("  Running %s...", method_name))
  
  result <- tryCatch({
    if (method_name == "xcell") {
      # Use xCell package directly (most reliable)
      scores <- as.data.frame(xCellAnalysis(exp_matrix))
      scores <- as.data.frame(t(scores))
      colnames(scores) <- paste0("xCell_", colnames(scores))
      scores
    } else if (method_name == "epic") {
      # Use EPIC from immunedeconv
      scores <- deconvolute(exp_matrix, "epic")
      scores_df <- as.data.frame(t(scores[, -1]))
      colnames(scores_df) <- paste0("EPIC_", scores$cell_type)
      scores_df
    } else if (method_name == "quantiseq") {
      # Use quanTIseq from immunedeconv
      scores <- deconvolute(exp_matrix, "quantiseq")
      scores_df <- as.data.frame(t(scores[, -1]))
      colnames(scores_df) <- paste0("quanTIseq_", scores$cell_type)
      scores_df
    } else if (method_name == "mcp_counter") {
      # Use MCPcounter from immunedeconv
      scores <- deconvolute(exp_matrix, "mcp_counter")
      scores_df <- as.data.frame(t(scores[, -1]))
      colnames(scores_df) <- paste0("MCPcounter_", scores$cell_type)
      scores_df
    } else if (method_name == "estimate") {
      # Use ESTIMATE to get immune/stromal scores
      library(estimate)
      # Prepare input file
      write.table(cbind(GeneSymbol=rownames(exp_matrix), exp_matrix),
                  file="temp_input.gct", sep="\t", quote=FALSE)
      
      # Run ESTIMATE
      filterCommonGenes(input.f="temp_input.gct", 
                        output.f="temp_filtered.gct",
                        id="GeneSymbol")
      estimateScore("temp_filtered.gct", "temp_estimate.gct")
      
      # Read results
      estimate_scores <- read.delim("temp_estimate.gct", skip=2, row.names=1)
      scores <- as.data.frame(t(estimate_scores[1:3, -1]))
      colnames(scores) <- c("ESTIMATE_Stromal", "ESTIMATE_Immune", "ESTIMATE_ESTIMATE")
      scores
    }
  }, error = function(e) {
    cat(sprintf(" FAILED: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(result)) {
    cat(" SUCCESS\n")
  }
  
  return(result)
}

# Run each method
methods <- c("xcell", "epic", "quantiseq", "mcp_counter", "estimate")
immune_results <- list()

for (method in methods) {
  immune_results[[method]] <- safe_deconvolute(method, exp_matrix)
}

# ============================================================================
# 5. COMBINE IMMUNE SCORES
# ============================================================================

cat("\nStep 3: Combining immune scores...\n")

# Combine all successful methods
immune_scores <- data.frame()

for (method_name in names(immune_results)) {
  if (!is.null(immune_results[[method_name]])) {
    scores_df <- immune_results[[method_name]]
    
    if (nrow(immune_scores) == 0) {
      immune_scores <- scores_df
    } else {
      common_samples <- intersect(rownames(immune_scores), rownames(scores_df))
      if (length(common_samples) > 0) {
        immune_scores <- cbind(immune_scores[common_samples, , drop = FALSE],
                               scores_df[common_samples, , drop = FALSE])
      }
    }
  }
}

# Check what we got
cat(sprintf("Successfully obtained immune scores from %d methods\n", 
            sum(!sapply(immune_results, is.null))))

if (nrow(immune_scores) == 0) {
  stop("No immune scores were obtained. Check your data format.")
}

# Select key immune cell types for breast cancer
cat("Available immune cell types:\n")
print(colnames(immune_scores))

# Prioritize these cell types (based on what's available)
target_cells <- c(
  # xCell cells
  "xCell_CD8+ T-cells", "xCell_CD4+ T-cells", "xCell_Tregs",
  "xCell_Macrophages", "xCell_Macrophages M1", "xCell_Macrophages M2",
  "xCell_B-cells", "xCell_NK cells", "xCell_Neutrophils",
  "xCell_Cancer associated fibroblasts", "xCell_Endothelial cells",
  
  # EPIC cells
  "EPIC_Bcells", "EPIC_CD4_Tcells", "EPIC_CD8_Tcells", "EPIC_Macrophages",
  "EPIC_NKcells", "EPIC_Endothelial", "EPIC_Fibroblasts",
  
  # quanTIseq cells
  "quanTIseq_B cells", "quanTIseq_CD8 T cells", "quanTIseq_Macrophages",
  
  # MCPcounter cells
  "MCPcounter_CD8 T cells", "MCPcounter_Monocytic lineage",
  "MCPcounter_Fibroblasts",
  
  # ESTIMATE scores
  "ESTIMATE_Immune", "ESTIMATE_Stromal"
)

# Filter to available cells
available_targets <- target_cells[target_cells %in% colnames(immune_scores)]
immune_scores <- immune_scores[, available_targets, drop = FALSE]

cat(sprintf("\nSelected %d immune cell types for analysis:\n", ncol(immune_scores)))
print(colnames(immune_scores))

# Save immune scores
write.csv(immune_scores, "immune_cell_scores.csv", row.names = TRUE)
cat("\nSaved immune scores to 'immune_cell_scores.csv'\n")

# ============================================================================
# 6. CORRELATION ANALYSIS
# ============================================================================

cat("\nStep 4: Calculating gene-immune cell correlations...\n")

# Prepare gene expression data
gene_exp_t <- as.data.frame(t(gene_exp))

# Ensure same samples in both datasets
common_samples <- intersect(rownames(gene_exp_t), rownames(immune_scores))
cat(sprintf("Using %d common samples for correlation analysis\n", length(common_samples)))

if (length(common_samples) < 10) {
  warning("Low sample count for correlation analysis. Results may not be reliable.")
}

gene_exp_t <- gene_exp_t[common_samples, , drop = FALSE]
immune_scores <- immune_scores[common_samples, , drop = FALSE]

# Calculate Spearman correlations
cor_results <- data.frame()

for (gene in colnames(gene_exp_t)) {
  cat(sprintf("  Calculating correlations for %s...\n", gene))
  
  for (cell_type in colnames(immune_scores)) {
    # Skip if too many missing values
    if (sum(!is.na(gene_exp_t[[gene]])) < 5 || 
        sum(!is.na(immune_scores[[cell_type]])) < 5) {
      next
    }
    
    # Calculate correlation
    cor_test <- tryCatch({
      cor.test(gene_exp_t[[gene]], 
               immune_scores[[cell_type]], 
               method = "spearman",
               exact = FALSE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(cor_test)) {
      cor_results <- rbind(cor_results, data.frame(
        Gene = gene,
        Immune_Cell = cell_type,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Adjust p-values for multiple testing
if (nrow(cor_results) > 0) {
  cor_results$FDR <- p.adjust(cor_results$P_value, method = "BH")
  
  # Sort by absolute correlation
  cor_results <- cor_results[order(-abs(cor_results$Correlation)), ]
  
  # Save results
  write.csv(cor_results, "gene_immune_correlations.csv", row.names = FALSE)
  cat(sprintf("Calculated %d correlations. Saved to file.\n", nrow(cor_results)))
} else {
  stop("No correlations could be calculated. Check your data.")
}

# ============================================================================
# 7. VISUALIZATION
# ============================================================================

cat("\nStep 5: Creating visualizations...\n")

# 7.1 Create output directory for plots
if (!dir.exists("plots11")) {
  dir.create("plots11")
}

# 7.2 Heatmap of top correlations
pdf("plots/correlation_heatmap.pdf", width = 12, height = 10)

# Select top correlations (up to 30)
top_n <- min(30, nrow(cor_results))
top_cors <- head(cor_results, top_n)

# Create matrix for heatmap
heatmap_data <- dcast(top_cors, Gene ~ Immune_Cell, value.var = "Correlation")
rownames(heatmap_data) <- heatmap_data$Gene
heatmap_matrix <- as.matrix(heatmap_data[, -1, drop = FALSE])

# Create heatmap
corrplot(heatmap_matrix,
         method = "color",
         type = "full",
         tl.col = "black",
         tl.cex = 0.7,
         tl.srt = 45,
         mar = c(0, 0, 2, 0),
         title = paste("Top", top_n, "Gene-Immune Cell Correlations"),
         addCoef.col = "black",
         number.cex = 0.6,
         cl.cex = 0.7,
         is.corr = FALSE,
         number.digits = 2)

dev.off()
cat("Created correlation_heatmap.pdf\n")

# 7.3 Scatter plots for top 6 correlations
pdf("plots/top_correlations_scatter.pdf", width = 14, height = 12)

# Get top 6 correlations
top_6 <- head(cor_results, 6)

# Create layout
layout_matrix <- matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
layout(layout_matrix)

for (i in 1:nrow(top_6)) {
  gene <- top_6$Gene[i]
  cell <- top_6$Immune_Cell[i]
  cor_val <- round(top_6$Correlation[i], 3)
  p_val <- formatC(top_6$P_value[i], format = "e", digits = 2)
  fdr_val <- formatC(top_6$FDR[i], format = "e", digits = 2)
  
  # Create scatter plot
  plot(gene_exp_t[[gene]], immune_scores[[cell]],
       xlab = paste(gene, "Expression"),
       ylab = paste(cell, "Score"),
       main = paste(gene, "vs", cell),
       pch = 19, col = rgb(0.2, 0.4, 0.8, 0.6),
       cex = 1.0)
  
  # Add regression line
  if (sum(!is.na(gene_exp_t[[gene]])) > 2 && 
      sum(!is.na(immune_scores[[cell]])) > 2) {
    lm_fit <- lm(immune_scores[[cell]] ~ gene_exp_t[[gene]])
    abline(lm_fit, col = "red", lwd = 2)
  }
  
  # Add correlation info
  mtext(paste0("r = ", cor_val, "\np = ", p_val, "\nFDR = ", fdr_val), 
        side = 3, line = -2, adj = 0.95, cex = 0.7)
}

dev.off()
cat("Created top_correlations_scatter.pdf\n")

# 7.4 Bar plot of strongest correlations per gene
pdf("plots/gene_immune_summary.pdf", width = 10, height = 8)

# Get strongest positive and negative correlation for each gene
gene_summary <- cor_results %>%
  group_by(Gene) %>%
  summarize(
    Max_Positive = max(Correlation, na.rm = TRUE),
    Max_Negative = min(Correlation, na.rm = TRUE),
    Pos_Cell = Immune_Cell[which.max(Correlation)],
    Neg_Cell = Immune_Cell[which.min(Correlation)],
    Num_Sig_Pos = sum(Correlation > 0.3 & FDR < 0.1, na.rm = TRUE),
    Num_Sig_Neg = sum(Correlation < -0.3 & FDR < 0.1, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(Max_Positive, Max_Negative),
               names_to = "Direction",
               values_to = "Correlation") %>%
  mutate(Cell_Type = ifelse(Direction == "Max_Positive", Pos_Cell, Neg_Cell),
         Cell_Short = gsub("^[^_]+_", "", Cell_Type))  # Remove method prefix

# Create bar plot
p <- ggplot(gene_summary, aes(x = Gene, y = Correlation, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Cell_Short, 
                y = Correlation + sign(Correlation) * 0.03),
            position = position_dodge(width = 0.8),
            size = 3, angle = 45, hjust = 0) +
  labs(title = "Strongest Immune Cell Correlations for Each Gene",
       subtitle = "Positive (red) and Negative (blue) Associations",
       y = "Correlation Coefficient (Spearman's rho)",
       x = "") +
  scale_fill_manual(values = c("Max_Positive" = "indianred2", 
                               "Max_Negative" = "steelblue2"),
                    labels = c("Most Positive", "Most Negative")) +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank())

# Add annotation for number of significant correlations
annotation_df <- gene_summary %>%
  group_by(Gene) %>%
  summarize(Num_Sig_Pos = first(Num_Sig_Pos),
            Num_Sig_Neg = first(Num_Sig_Neg))

p <- p + geom_text(data = annotation_df, 
                   aes(x = Gene, y = -0.9, 
                       label = paste0("+", Num_Sig_Pos, "/-", Num_Sig_Neg)),
                   inherit.aes = FALSE, size = 3, color = "gray40")

print(p)
dev.off()
cat("Created gene_immune_summary.pdf\n")

# 7.5 Dot plot of all significant correlations
pdf("plots/significant_correlations_dotplot.pdf", width = 12, height = 8)

# Filter significant correlations
sig_cors <- cor_results %>%
  filter(FDR < 0.1) %>%
  arrange(Correlation)

if (nrow(sig_cors) > 0) {
  # Create dot plot
  p <- ggplot(sig_cors, aes(x = Correlation, y = reorder(Immune_Cell, Correlation), 
                            color = Gene, size = abs(Correlation))) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Significant Gene-Immune Cell Correlations (FDR < 0.1)",
         x = "Correlation Coefficient",
         y = "Immune Cell Type") +
    scale_color_brewer(palette = "Set1") +
    scale_size_continuous(range = c(2, 6)) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          legend.position = "right")
  
  print(p)
} else {
  plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "No significant correlations (FDR < 0.1)", cex = 1.2)
}

dev.off()
cat("Created significant_correlations_dotplot.pdf\n")

# ============================================================================
# 8. GENERATE COMPREHENSIVE REPORT
# ============================================================================

cat("\nStep 6: Generating comprehensive report...\n")

# Create summary statistics
summary_stats <- cor_results %>%
  group_by(Gene) %>%
  summarize(
    Num_Total_Correlations = n(),
    Num_Sig_Positive = sum(Correlation > 0 & FDR < 0.1, na.rm = TRUE),
    Num_Sig_Negative = sum(Correlation < 0 & FDR < 0.1, na.rm = TRUE),
    Mean_Correlation = mean(Correlation, na.rm = TRUE),
    Median_Correlation = median(Correlation, na.rm = TRUE),
    Strongest_Positive = max(Correlation, na.rm = TRUE),
    Strongest_Negative = min(Correlation, na.rm = TRUE),
    Most_Positive_Cell = Immune_Cell[which.max(Correlation)],
    Most_Negative_Cell = Immune_Cell[which.min(Correlation)],
    Top_Immunosuppressive = ifelse(
      any(grepl("M2|Treg|Fibroblast", Immune_Cell) & FDR < 0.1 & Correlation > 0.3),
      paste(Immune_Cell[grepl("M2|Treg|Fibroblast", Immune_Cell) & 
                          FDR < 0.1 & 
                          Correlation > 0.3], collapse = "; "),
      "None"
    )
  )

# Save summary
write.csv(summary_stats, "correlation_summary_statistics.csv", row.names = FALSE)

# Generate HTML report
sink("analysis_report.html")

cat('<!DOCTYPE html>
<html>
<head>
<title>Immune Reprogramming Analysis Report</title>
<style>
body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
h1, h2, h3 { color: #2c3e50; }
table { border-collapse: collapse; width: 100%; margin: 20px 0; }
th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
th { background-color: #f2f2f2; }
tr:nth-child(even) { background-color: #f9f9f9; }
.section { margin: 30px 0; padding: 20px; background-color: #f8f9fa; border-radius: 5px; }
.highlight { background-color: #fff3cd; padding: 10px; border-radius: 3px; }
</style>
</head>
<body>')

cat('<h1>Immune Reprogramming Analysis Report</h1>')
cat('<p><strong>Date:</strong>', format(Sys.Date(), '%Y-%m-%d'), '</p>')
cat('<p><strong>Time:</strong>', format(Sys.time(), '%H:%M:%S'), '</p>')

cat('<div class="section">
<h2>Analysis Summary</h2>
<p><strong>Genes analyzed:</strong>', paste(available_genes, collapse=', '), '</p>
<p><strong>Samples analyzed:</strong>', nrow(immune_scores), '</p>
<p><strong>Immune cell types evaluated:</strong>', ncol(immune_scores), '</p>
</div>')

cat('<div class="section">
<h2>Key Findings</h2>')

# Top 5 correlations
top_5 <- head(cor_results, 5)
cat('<h3>Top 5 Strongest Correlations</h3>
<table>
<tr><th>Gene</th><th>Immune Cell</th><th>Correlation</th><th>P-value</th><th>FDR</th></tr>')

for(i in 1:nrow(top_5)) {
  cat('<tr>')
  cat('<td>', top_5$Gene[i], '</td>')
  cat('<td>', top_5$Immune_Cell[i], '</td>')
  cat('<td>', round(top_5$Correlation[i], 3), '</td>')
  cat('<td>', formatC(top_5$P_value[i], format='e', digits=2), '</td>')
  cat('<td>', formatC(top_5$FDR[i], format='e', digits=2), '</td>')
  cat('</tr>')
}
cat('</table>')

cat('<div class="highlight">
<h3>Interpretation Guide</h3>
<ul>
<li><strong>Positive correlation with M2 macrophages/Tregs/CAFs</strong>: Suggests immunosuppressive role</li>
<li><strong>Negative correlation with CD8+ T cells/NK cells</strong>: Suggests immune evasion role</li>
<li><strong>Correlation > |0.4| with FDR < 0.1</strong>: Strong evidence of association</li>
</ul>
</div>')

cat('</div>')

cat('<div class="section">
<h2>Gene-Specific Results</h2>')

for(gene in available_genes) {
  gene_data <- cor_results[cor_results$Gene == gene, ]
  sig_pos <- gene_data[gene_data$Correlation > 0.3 & gene_data$FDR < 0.1, ]
  sig_neg <- gene_data[gene_data$Correlation < -0.3 & gene_data$FDR < 0.1, ]
  
  cat('<h3>', gene, '</h3>')
  cat('<p><strong>Total correlations calculated:</strong>', nrow(gene_data), '</p>')
  
  if(nrow(sig_pos) > 0) {
    cat('<p><strong>Significant positive associations (r > 0.3, FDR < 0.1):</strong></p><ul>')
    for(i in 1:min(3, nrow(sig_pos))) {
      cat('<li>', sig_pos$Immune_Cell[i], ' (r = ', round(sig_pos$Correlation[i], 2), 
          ', FDR = ', formatC(sig_pos$FDR[i], format='e', digits=1), ')</li>', sep='')
    }
    cat('</ul>')
  }
  
  if(nrow(sig_neg) > 0) {
    cat('<p><strong>Significant negative associations (r < -0.3, FDR < 0.1):</strong></p><ul>')
    for(i in 1:min(3, nrow(sig_neg))) {
      cat('<li>', sig_neg$Immune_Cell[i], ' (r = ', round(sig_neg$Correlation[i], 2), 
          ', FDR = ', formatC(sig_neg$FDR[i], format='e', digits=1), ')</li>', sep='')
    }
    cat('</ul>')
  }
}

cat('</div>')

cat('<div class="section">
<h2>Generated Files</h2>
<ul>
<li><strong>immune_cell_scores.csv</strong>: Immune cell abundance scores</li>
<li><strong>gene_immune_correlations.csv</strong>: All correlation results</li>
<li><strong>correlation_summary_statistics.csv</strong>: Summary statistics</li>
<li><strong>plots/correlation_heatmap.pdf</strong>: Heatmap visualization</li>
<li><strong>plots/top_correlations_scatter.pdf</strong>: Scatter plots</li>
<li><strong>plots/gene_immune_summary.pdf</strong>: Summary bar plot</li>
<li><strong>plots/significant_correlations_dotplot.pdf</strong>: Dot plot of significant correlations</li>
</ul>
</div>')

cat('<div class="highlight">
<h2>Next Steps</h2>
<ol>
<li>Examine the generated plots to visualize correlations</li>
<li>Focus on genes showing strong correlations (|r| > 0.4) with immune cells</li>
<li>Prioritize associations with immunosuppressive cells (M2 macrophages, Tregs, CAFs)</li>
<li>Consider experimental validation for top findings</li>
<li>For deeper analysis, consider single-cell RNA-seq validation</li>
</ol>
</div>')

cat('</body></html>')

sink()

# Also print summary to console
cat("\n" + strrep("=", 70) + "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 70) + "\n\n")

cat("SUCCESSFULLY GENERATED:\n")
cat("1. Immune cell scores: immune_cell_scores.csv\n")
cat("2. Correlation results: gene_immune_correlations.csv\n")
cat("3. Summary statistics: correlation_summary_statistics.csv\n")
cat("4. Visualizations in 'plots/' directory:\n")
cat("   - correlation_heatmap.pdf\n")
cat("   - top_correlations_scatter.pdf\n")
cat("   - gene_immune_summary.pdf\n")
cat("   - significant_correlations_dotplot.pdf\n")
cat("5. HTML report: analysis_report.html\n\n")

cat("KEY STATISTICS:\n")
cat(sprintf("- Samples analyzed: %d\n", nrow(immune_scores)))
cat(sprintf("- Immune cell types: %d\n", ncol(immune_scores)))
cat(sprintf("- Total correlations: %d\n", nrow(cor_results)))
cat(sprintf("- Significant correlations (FDR < 0.1): %d\n", 
            sum(cor_results$FDR < 0.1, na.rm = TRUE)))
cat(sprintf("- Strong correlations (|r| > 0.4): %d\n", 
            sum(abs(cor_results$Correlation) > 0.4, na.rm = TRUE)))

cat("\nTOP 3 FINDINGS:\n")
for(i in 1:min(3, nrow(cor_results))) {
  row <- cor_results[i, ]
  cat(sprintf("%d. %s vs %s: r = %.2f (FDR = %.2e)\n", 
              i, row$Gene, row$Immune_Cell, row$Correlation, row$FDR))
}

cat("\n" + strrep("=", 70) + "\n")