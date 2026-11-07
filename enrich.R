#1. Setup and gene list preparation
############################################################
# ASSUMED INPUT
#   deg: data.frame from limma::topTable (all genes)
#        with columns: logFC, adj.P.Val, etc.
############################################################

# If you already have 'deg' in memory from the limma pipeline, skip this.
# Otherwise, for example:
deg <- read.csv("DEG_full_results.csv", row.names = 1)

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define significant DEGs (you can tweak thresholds)
deg_sig <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Extract gene symbols
gene_symbols <- rownames(deg_sig)

# Map to Entrez IDs
gene_map <- bitr(gene_symbols,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

entrez_deg <- unique(gene_map$ENTREZID)
##
#2. GO Biological Process enrichment
##
ego_bp <- enrichGO(
  gene          = entrez_deg,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# View top GO terms
head(as.data.frame(ego_bp), 20)

# Save
write.csv(as.data.frame(ego_bp), "GO_BP_enrichment_DEGs.csv", row.names = FALSE)
##
#3. KEGG pathway enrichment
##
ekegg <- enrichKEGG(
  gene          = entrez_deg,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Make gene IDs human-readable
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# View top KEGG pathways
head(as.data.frame(ekegg), 20)

# Save
write.csv(as.data.frame(ekegg), "KEGG_enrichment_DEGs.csv", row.names = FALSE)

##=============================================================
# gene selection
##=============================================================

# deg_sig: significant DEGs from limma
# ego_bp: GO enrichment result
# ekegg: KEGG enrichment result

library(dplyr)

# 1. Start with significant DEGs
deg_sig <- deg %>% filter(adj.P.Val < 0.001 & abs(logFC) > 1.2)

# 2. Extract pathway genes from GO and KEGG
go_genes <- unique(unlist(strsplit(ego_bp$geneID, "/")))
kegg_genes <- unique(unlist(strsplit(ekegg$geneID, "/")))

# 3. Keep DEGs that are in enriched pathways
candidate_genes <- intersect(rownames(deg_sig), union(go_genes, kegg_genes))

# 4. Remove redundant genes by correlation filtering
expr_sub <- expr[candidate_genes, ]
cor_mat <- cor(t(log2(expr_sub + 1)))

# Keep one gene per highly correlated cluster (|r| > 0.9)
keep <- c()
used <- c()

for (g in rownames(cor_mat)) {
  if (!(g %in% used)) {
    cluster <- names(which(abs(cor_mat[g, ]) > 0.9))
    keep <- c(keep, g)
    used <- c(used, cluster)
  }
}

selected_genes <- keep
selected_genes
##===============================================================
# Refine selected genes
##===============================================================
log_expr <- log2(expr[selected_genes, ] + 1)
cor_mat <- cor(t(log_expr))
dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "average")
clusters <- cutree(hc, k = 25)
rep_genes <- sapply(unique(clusters), function(cl){
  genes <- names(clusters[clusters == cl])
  sub <- deg[genes, ]
  genes[which.min(sub$adj.P.Val)]
})
rep_genes <- unname(rep_genes)
rep_genes1=c(rep_genes, "COL3A1", "FN1", "LOX", "ITGB1", "ITGAV", "ARHGAP5", "DOCK9", "TGFB2", "TGFBR2", "RABEP1", "ARFGEF1")
rep_genes=rep_genes1

##Elastic net
library(glmnet)

x <- t(log2(expr[rep_genes, ] + 1))
y <- factor(group)

fit <- cv.glmnet(x, y, family = "binomial", alpha = 0.5)
coef_genes <- rownames(coef(fit, s = "lambda.min"))[coef(fit, s = "lambda.min") != 0]
# Extract coefficient object
cf <- coef(fit, s = "lambda.min")

# Convert to a numeric vector
cf_vec <- as.numeric(cf)

# Extract gene names (rownames of the sparse matrix)
cf_names <- rownames(cf)

# Keep only non-zero coefficients
selected_genes_glmnet <- cf_names[cf_vec != 0]
selected_genes_glmnet <- selected_genes_glmnet[selected_genes_glmnet != "(Intercept)"]
coef_df <- data.frame(
  gene = cf_names,
  coef = cf_vec
)

coef_df <- coef_df[coef_df$gene != "(Intercept)" & coef_df$coef != 0, ]
coef_df <- coef_df[order(abs(coef_df$coef), decreasing = TRUE), ]

coef_df

write.csv(coef_df,file="coef_df.csv",row.names = F)
