# Required packages
install_if_missing <- function(pkgs){
  to_install <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
}
install_if_missing(c("tidyverse","limma","pheatmap","reshape2","ggpubr","data.table"))

library(tidyverse)
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(data.table)

# -------------------------
# User inputs: update paths
# -------------------------
setwd("D:/GSVAMODEL")
gsva_file      <- "gsvaout.txt"        # rows = pathway id, cols = samples; first col named "pathway" or "id"
expr_file      <- "TCGA_TPM.txt"          # rows = gene symbol, cols = samples; first col named "gene" or "symbol"
meta_file      <- "C"           # must contain columns: sample, subtype (values: TypeA, TypeC)
gmt_file       <- "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt"           # optional: GMT file with gene sets (if available)
out_prefix     <- "gsva_top6"              # prefix for saved plots/tables (optional)

# Which contrast: set test and reference labels exactly as in metadata
test_group     <- "3"
ref_group      <- "1"

# Keywords to select cell-cycle / condensin related pathways if you want automatic selection
pathway_keywords <- c("CELL_CYCLE","CONDENSIN","DNA_REPLICATION","G1_S","G2_M","P27","P16","CDC25","WEE1","REPLICATION")

# -------------------------
# Read data
# -------------------------
gsva <- fread(gsva_file, data.table = FALSE)
# ensure first column is pathway id
if(!("pathway" %in% colnames(gsva)) & !("id" %in% colnames(gsva))){
  colnames(gsva)[1] <- "pathway"
}
rownames(gsva) <- gsva$pathway
gsva <- gsva %>% select(-pathway)

expr <- fread(expr_file, data.table = FALSE)
if(!("gene" %in% colnames(expr)) & !("symbol" %in% colnames(expr))){
  colnames(expr)[1] <- "gene"
}
rownames(expr) <- expr[[1]]
expr <- expr %>% select(-1)

meta <- fread(meta_file, data.table = FALSE)
if(!("sample" %in% colnames(meta))) stop("metadata must contain a 'sample' column matching sample IDs")
if(!("subtype" %in% colnames(meta))) stop("metadata must contain a 'subtype' column with TypeA/TypeC labels")

# Ensure sample order matches across matrices
samples_common <- Reduce(intersect, list(colnames(gsva), colnames(expr), meta$sample))
if(length(samples_common) < 3) stop("Too few matching samples across files. Check sample IDs.")
gsva <- gsva[, samples_common, drop = FALSE]
expr <- expr[, samples_common, drop = FALSE]
meta <- meta %>% filter(sample %in% samples_common) %>% arrange(match(sample, samples_common))

# -------------------------
# If you have a GMT file, read it; otherwise user must supply gene lists for pathways
# -------------------------
read_gmt <- function(gmt){
  if(!file.exists(gmt)) return(NULL)
  lines <- readLines(gmt)
  gsets <- lapply(lines, function(l){
    parts <- strsplit(l, "\t")[[1]]
    name <- parts[1]
    genes <- parts[-c(1,2)]
    list(name = name, genes = genes)
  })
  names(gsets) <- sapply(gsets, `[[`, "name")
  return(gsets)
}
gsets <- read_gmt(gmt_file)

# -------------------------
# Select top 6 cell-cycle/condensin pathways
# Option A: if you already have a limma results table (adj.P.Val etc), read it and pick top6
# Option B: pick pathways by keyword and then by variance or mean difference
# -------------------------
# Try to auto-select by keyword first
pw_names <- rownames(gsva)
sel_idx <- unique(unlist(lapply(pathway_keywords, function(k) grep(k, pw_names, ignore.case = TRUE))))
if(length(sel_idx) >= 6){
  candidate_pws <- pw_names[sel_idx]
  # rank by variance across samples (as proxy) and pick top 6
  pw_var <- apply(gsva[candidate_pws, , drop = FALSE], 1, var)
  top6 <- names(sort(pw_var, decreasing = TRUE))[1:6]
} else {
  # fallback: pick top 6 by overall variance across all pathways
  pw_var_all <- apply(gsva, 1, var)
  top6 <- names(sort(pw_var_all, decreasing = TRUE))[1:6]
}
top6

# -------------------------
# Plot per-sample GSVA scores for the top6
# -------------------------
gsva_long <- gsva[top6, ] %>% 
  as.data.frame() %>% 
  rownames_to_column("pathway") %>% 
  pivot_longer(-pathway, names_to = "sample", values_to = "gsva")

gsva_long <- gsva_long %>% left_join(meta %>% rename(sample = sample), by = "sample")
gsva_long$pathway <- factor(gsva_long$pathway, levels = top6)

p_list <- lapply(top6, function(pw){
  df <- gsva_long %>% filter(pathway == pw)
  ggplot(df, aes(x = subtype, y = gsva, color = subtype)) +
    geom_violin(trim = FALSE, alpha = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
    stat_summary(fun = median, geom = "point", size = 2, color = "black") +
    theme_minimal() +
    labs(title = pw, x = "", y = "GSVA score")
})

# arrange and save
plot_grid <- ggarrange(plotlist = p_list, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave(paste0(out_prefix, "_top6_gsva_violin.png"), plot_grid, width = 14, height = 8, dpi = 300)

# -------------------------
# Recompute limma on GSVA scores for these pathways (TypeA - TypeC)
# -------------------------
design <- model.matrix(~ 0 + factor(meta$subtype))
colnames(design) <- gsub("factor\\(meta\\$subtype\\)", "", colnames(design))
if(!(test_group %in% colnames(design)) | !(ref_group %in% colnames(design))){
  stop("Test or reference group not found in metadata subtype column. Check labels.")
}
contrast <- makeContrasts(contr = paste0(test_group, "-", ref_group), levels = design)
fit <- lmFit(as.matrix(gsva), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
gsva_res <- topTable(fit2, number = nrow(gsva), sort.by = "P")
write.csv(gsva_res, paste0(out_prefix, "_gsva_limma_results.csv"), row.names = TRUE)

# Subset results for top6
gsva_res_top6 <- gsva_res[rownames(gsva_res) %in% top6, ]
gsva_res_top6 <- gsva_res_top6[match(top6, rownames(gsva_res_top6)), ]
gsva_res_top6

# -------------------------
# For each pathway: top 5 driving genes
# Two methods:
# 1) Correlation of each gene's expression with pathway GSVA score (Pearson)
# 2) Leading-edge proxy: among pathway member genes, rank by differential expression (limma on gene expression)
# -------------------------
# 1) Correlation method
corr_top5 <- list()
for(pw in top6){
  pw_scores <- as.numeric(gsva[pw, ])
  cors <- apply(expr, 1, function(g) cor(g, pw_scores, use = "pairwise.complete.obs"))
  cors_df <- data.frame(gene = rownames(expr), cor = cors, stringsAsFactors = FALSE)
  cors_df <- cors_df %>% arrange(desc(abs(cor)))
  corr_top5[[pw]] <- head(cors_df, 5)
}

# 2) Leading-edge proxy: differential expression (TypeA - TypeC) on genes
# Build design for gene expression
design_g <- model.matrix(~ 0 + factor(meta$subtype))
colnames(design_g) <- gsub("factor\\(meta\\$subtype\\)", "", colnames(design_g))
fit_g <- lmFit(as.matrix(expr), design_g)
fit_g2 <- contrasts.fit(fit_g, contrast)
fit_g2 <- eBayes(fit_g2)
gene_res <- topTable(fit_g2, number = nrow(expr), sort.by = "P")
rownames(gene_res) <- rownames(expr)

lead_top5 <- list()
for(pw in top6){
  # get gene list for pathway from GMT if available, else try to parse pathway name to gene list (not possible)
  if(!is.null(gsets) && pw %in% names(gsets)){
    members <- gsets[[pw]]$genes
  } else {
    # fallback: try to find genes by splitting pathway name (not reliable)
    members <- intersect(rownames(expr), character(0))
  }
  if(length(members) == 0){
    lead_top5[[pw]] <- data.frame(gene = NA, logFC = NA, t = NA, P.Value = NA)[0,]
  } else {
    members_in_expr <- intersect(members, rownames(gene_res))
    if(length(members_in_expr) == 0){
      lead_top5[[pw]] <- data.frame(gene = NA, logFC = NA, t = NA, P.Value = NA)[0,]
    } else {
      sub <- gene_res[members_in_expr, , drop = FALSE]
      sub <- sub %>% rownames_to_column("gene") %>% arrange(P.Value)
      lead_top5[[pw]] <- head(sub %>% select(gene, logFC, t, P.Value), 5)
    }
  }
}

# Save results
saveRDS(list(corr_top5 = corr_top5, lead_top5 = lead_top5), file = paste0(out_prefix, "_top5_genes_lists.rds"))
# Also write human-readable CSVs
for(pw in top6){
  write.csv(corr_top5[[pw]], paste0(out_prefix, "_", make.names(pw), "_corr_top5.csv"), row.names = FALSE)
  if(nrow(lead_top5[[pw]])>0) write.csv(lead_top5[[pw]], paste0(out_prefix, "_", make.names(pw), "_lead_top5.csv"), row.names = FALSE)
}

# -------------------------
# Heatmap of union of top genes across the 6 sets (correlation top genes)
# -------------------------
union_genes <- unique(unlist(lapply(corr_top5, function(x) x$gene)))
expr_sub <- expr[intersect(union_genes, rownames(expr)), , drop = FALSE]
# scale rows for visualization
expr_scaled <- t(scale(t(expr_sub)))
annotation_col <- data.frame(subtype = meta$subtype)
rownames(annotation_col) <- meta$sample
pheatmap(expr_scaled, annotation_col = annotation_col, show_rownames = TRUE,
         filename = paste0(out_prefix, "_top_genes_heatmap.png"), width = 8, height = 6)

# -------------------------
# Short diagnostics printed to console
# -------------------------
cat("Top6 pathways selected:\n")
print(top6)
cat("\nGSVA limma results for top6 (TypeA - TypeC):\n")
print(gsva_res_top6)
cat("\nTop 5 correlated genes per pathway (correlation method):\n")
print(lapply(corr_top5, head))
cat("\nTop 5 leading-edge proxy genes per pathway (if GMT provided):\n")
print(lead_top5)