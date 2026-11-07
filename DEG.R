############################################################
# FULL DEG ANALYSIS PIPELINE FOR TPM DATA (READY TO RUN)
############################################################

# Required packages
library(edgeR)
library(limma)
library(ggplot2)

############################################################
# 1. INPUTS
############################################################

expr=read.table("TCGA_TPM.txt",row.names = 1, check.names = F, sep = "\t",header=T )
label=read.csv("group_labels.csv")
group <-lebel$group
expr=expr[label$sample]

group <- factor(group)

############################################################
# 2. PREPROCESSING
############################################################

# Log-transform TPM
logTPM <- log2(expr + 0.1)

# Remove zero-variance genes
gene_var <- apply(logTPM, 1, var)
logTPM_filt <- logTPM[gene_var > 0, ]

# Convert logTPM back to pseudo-counts for voom
pseudo_counts <- 2^logTPM_filt

############################################################
# 3. CREATE DGEList AND DESIGN MATRIX
############################################################

dge <- DGEList(counts = pseudo_counts)
design <- model.matrix(~ group)

############################################################
# 4. VOOM TRANSFORMATION
############################################################

v <- voom(dge, design, plot = TRUE)

############################################################
# 5. LINEAR MODEL FITTING
############################################################

fit <- lmFit(v, design)
fit <- eBayes(fit)

############################################################
# 6. EXTRACT DEGs
############################################################

deg <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# Save significant DEGs
deg_sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]

write.csv(deg, "DEG_full_results.csv")
write.csv(deg_sig, "DEG_significant.csv")

############################################################
# 7. VOLCANO PLOT
############################################################

deg$gene <- rownames(deg)

ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )
