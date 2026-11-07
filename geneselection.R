############################################################
# 0. INPUTS
############################################################
# expr: gene x sample expression matrix (TPM or logTPM)
# group: factor of sample labels (two groups)
# deg: full limma DEG table
# ego_bp: GO BP enrichment result
# ekegg: KEGG enrichment result
expr=read.table("TCGA_TPM.txt",row.names = 1, check.names = F, sep = "\t",header=T )
label=read.csv("group_labels.csv")
group <-label$group
expr=expr[label$sample]
group <- factor(group)
deg=read.table("DEG_full_results.csv", row.names = 1,check.names = F, sep = ",",header=T )
ego_bp=read.table("GO_BP_enrichment_DEGs.csv", row.names = 1,check.names = F, sep = ",",header=T )
ekegg=read.table("KEGG_enrichment_DEGs.csv", row.names = 3,check.names = F, sep = ",",header=T )
library(dplyr)
library(glmnet)

############################################################
# 1. SELECT BIOLOGICALLY RELEVANT GENES (DEG ∩ GO/KEGG)
############################################################

# Significant DEGs
deg_sig <- deg %>% filter(adj.P.Val < 0.01)

# Extract GO/KEGG gene lists
go_genes <- unique(unlist(strsplit(ego_bp$geneID, "/")))
kegg_genes <- unique(unlist(strsplit(ekegg$geneID, "/")))

# Candidate genes = DEGs that appear in enriched pathways
candidate_genes <- intersect(rownames(deg_sig), union(go_genes, kegg_genes))

# Subset expression matrix
expr_sub <- expr[candidate_genes, ]

############################################################
# 2. CONVERT EXPRESSION TO RANKS (WITHIN EACH SAMPLE)
############################################################

rank_expr <- apply(expr_sub, 2, rank, ties.method = "average")

# rank_expr is now genes x samples (same shape as expr_sub)

############################################################
# 3. COMPUTE RANK-BASED DISCRIMINATIVE SCORE
############################################################
# For each gene, compute:
#   score = |mean(rank in group1) - mean(rank in group2)|

group1 <- which(group == levels(group)[1])
group2 <- which(group == levels(group)[2])

rank_score <- apply(rank_expr, 1, function(x){
  abs(mean(x[group1]) - mean(x[group2]))
})

# Keep top 200 rank-discriminative genes
top_rank_genes <- names(sort(rank_score, decreasing = TRUE))[1:200]

rank_expr_top <- rank_expr[top_rank_genes, ]

############################################################
# 4. FIT RANK-BASED ELASTIC NET CLASSIFIER
############################################################

x <- t(rank_expr_top)   # samples x genes
y <- group

fit <- cv.glmnet(
  x, y,
  family = "binomial",
  alpha = 0.5,        # elastic net
  standardize = FALSE # ranks should NOT be standardized
)

# Extract selected genes
cf <- coef(fit, s = "lambda.min")
cf_vec <- as.numeric(cf)
cf_names <- rownames(cf)

selected_genes_rank <- cf_names[cf_vec != 0]
selected_genes_rank <- selected_genes_rank[selected_genes_rank != "(Intercept)"]

selected_genes_rank
write.csv(selected_genes_rank, "selected_genes_rank.csv")
############################################################
# selected_genes_rank is your final rank-based signature
############################################################
