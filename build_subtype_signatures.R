# build_subtype_signatures.R
# Requirements: data.table, tidyverse, limma
# Run: Rscript build_subtype_signatures.R
setwd("D:/GSVAMODEL")
# ---------------------------
# 0. Parameters - edit paths
# ---------------------------
expr_file   <- "expr_log_tpm.csv"                       # genes x samples, first col = gene
gsva_file   <- "gsva_scores.csv"                        # pathways x samples, first col = pathway
meta_file   <- "meta.csv"                               # columns: sample, subtype (TypeA/TypeB/TypeC)
gmt_file    <- "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt"
diff_BA     <- "B-A_diff.csv"                           # B vs A (B - A)
diff_CA     <- "C-A_diff.csv"                           # C vs A (C - A)
diff_CB     <- "C-B_diff.csv"                           # C vs B (C - B)
out_prefix  <- "subtype_signatures"                     # output prefix
top_genes_per_pathway <- 10                             # how many leading genes to extract per pathway
signature_size <- 50                                    # final genes per subtype (adjustable)
adjp_thresh <- 0.05                                     # significance threshold for pathway diffs (optional filter)

# ---------------------------
# 1. Load libraries
# ---------------------------
install_if_missing <- function(pkgs){
  to_install <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
}
install_if_missing(c("data.table","tidyverse","limma"))
library(data.table)
library(tidyverse)
library(limma)

# ---------------------------
# 2. Read inputs and basic checks
# ---------------------------
expr <- fread(expr_file, data.table = FALSE)
if(!("gene" %in% colnames(expr))) colnames(expr)[1] <- "gene"
rownames(expr) <- expr$gene
expr <- expr[ , -1, drop = FALSE]

gsva <- fread(gsva_file, data.table = FALSE)
if(!("pathway" %in% colnames(gsva))) colnames(gsva)[1] <- "pathway"
rownames(gsva) <- gsva$pathway
gsva <- gsva[ , -1, drop = FALSE]

meta <- fread(meta_file, data.table = FALSE)
if(!all(c("sample","subtype") %in% colnames(meta))) stop("meta.csv must contain columns 'sample' and 'subtype'")

# sample intersection and ordering
common_samples <- Reduce(intersect, list(colnames(expr), colnames(gsva), meta$sample))
if(length(common_samples) < 3) stop("Too few matching samples across files. Check sample IDs.")
expr <- expr[, common_samples, drop = FALSE]
gsva <- gsva[, common_samples, drop = FALSE]
meta <- meta %>% filter(sample %in% common_samples) %>% arrange(match(sample, common_samples))
meta$subtype <- as.factor(meta$subtype)

# read diff tables
read_diff <- function(path){
  df <- fread(path, data.table = FALSE)
  if(!("id" %in% colnames(df))) colnames(df)[1] <- "id"
  rownames(df) <- df$id
  return(df)
}
diff_BA_df <- read_diff(diff_BA)   # B - A
diff_CA_df <- read_diff(diff_CA)   # C - A
diff_CB_df <- read_diff(diff_CB)   # C - B

# read GMT
read_gmt <- function(gmt){
  if(!file.exists(gmt)) return(NULL)
  lines <- readLines(gmt)
  gsets <- list()
  for(l in lines){
    parts <- strsplit(l, "\t")[[1]]
    name <- parts[1]
    genes <- parts[-c(1,2)]
    gsets[[name]] <- unique(genes)
  }
  return(gsets)
}
gsets <- read_gmt(gmt_file)
if(is.null(gsets)) message("No GMT found or could not read. DE-based leading edge will be limited.")

# ---------------------------
# 3. Define subtype-specific pathway lists
# Logic using provided diffs:
# - A > B and A > C  <=>  (B - A) < 0  AND  (C - A) < 0
# - B > A and B > C  <=>  (B - A) > 0  AND  (C - B) < 0
# - C > A and C > B  <=>  (C - A) > 0  AND  (C - B) > 0
# Optionally filter by adj.P.Val threshold in the diff tables.
# ---------------------------
pw_all <- rownames(gsva)

# helper to get numeric logFC from diff table (if missing, NA)
get_logfc <- function(df, pw) {
  if(is.null(df) || !(pw %in% rownames(df))) return(NA_real_)
  return(as.numeric(df[pw, "logFC"]))
}
get_adj <- function(df, pw) {
  if(is.null(df) || !(pw %in% rownames(df))) return(NA_real_)
  if("adj.P.Val" %in% colnames(df)) return(as.numeric(df[pw, "adj.P.Val"]))
  if("adj.P.Val" %in% colnames(df)) return(as.numeric(df[pw, "adj.P.Val"]))
  if("adj.P.Val" %in% colnames(df)) return(as.numeric(df[pw, "adj.P.Val"]))
  return(NA_real_)
}

pw_A <- pw_B <- pw_C <- character(0)
for(pw in pw_all){
  lfc_BA <- get_logfc(diff_BA_df, pw)   # B - A
  lfc_CA <- get_logfc(diff_CA_df, pw)   # C - A
  lfc_CB <- get_logfc(diff_CB_df, pw)   # C - B
  
  adj_BA <- get_adj(diff_BA_df, pw)
  adj_CA <- get_adj(diff_CA_df, pw)
  adj_CB <- get_adj(diff_CB_df, pw)
  
  # optional significance filter: require adj < adjp_thresh in both comparisons used
  sig_A <- (!is.na(adj_BA) && adj_BA < adjp_thresh) && (!is.na(adj_CA) && adj_CA < adjp_thresh)
  sig_B <- (!is.na(adj_BA) && adj_BA < adjp_thresh) && (!is.na(adj_CB) && adj_CB < adjp_thresh)
  sig_C <- (!is.na(adj_CA) && adj_CA < adjp_thresh) && (!is.na(adj_CB) && adj_CB < adjp_thresh)
  
  # Decide membership; if adjp present use significance, otherwise use logFC direction only
  if(!is.na(lfc_BA) && !is.na(lfc_CA)){
    if(lfc_BA < 0 && lfc_CA < 0 && (is.na(adj_BA) || is.na(adj_CA) || sig_A)) pw_A <- c(pw_A, pw)
    if(lfc_BA > 0 && !is.na(lfc_CB) && lfc_CB < 0 && (is.na(adj_BA) || is.na(adj_CB) || sig_B)) pw_B <- c(pw_B, pw)
    if(lfc_CA > 0 && !is.na(lfc_CB) && lfc_CB > 0 && (is.na(adj_CA) || is.na(adj_CB) || sig_C)) pw_C <- c(pw_C, pw)
  }
}

pw_A <- unique(pw_A); pw_B <- unique(pw_B); pw_C <- unique(pw_C)

# Save pathway lists
write.csv(data.frame(pathway = pw_A), paste0(out_prefix, "_pathways_TypeA.csv"), row.names = FALSE)
write.csv(data.frame(pathway = pw_B), paste0(out_prefix, "_pathways_TypeB.csv"), row.names = FALSE)
write.csv(data.frame(pathway = pw_C), paste0(out_prefix, "_pathways_TypeC.csv"), row.names = FALSE)

cat("Pathway counts: TypeA =", length(pw_A), "TypeB =", length(pw_B), "TypeC =", length(pw_C), "\n")

# ---------------------------
# 4. Leading-edge extraction functions
# ---------------------------
# 4A. Correlation-based leading edge (top N by absolute Pearson r)
leading_edge_cor <- function(pathway_name, gsva_mat, expr_mat, top_n = 30){
  if(!(pathway_name %in% rownames(gsva_mat))) return(NULL)
  scores <- as.numeric(gsva_mat[pathway_name, ])
  cors <- apply(expr_mat, 1, function(g) cor(g, scores, use = "pairwise.complete.obs"))
  df <- data.frame(gene = rownames(expr_mat), cor = cors, stringsAsFactors = FALSE)
  df <- df %>% arrange(desc(abs(cor))) %>% head(top_n)
  return(df)
}

# 4B. DE-based leading edge (pathway members only) using limma
leading_edge_de <- function(pathway_genes, expr_mat, meta_df, contrast_vec, top_n = 30){
  members <- intersect(pathway_genes, rownames(expr_mat))
  if(length(members) == 0) return(NULL)
  design <- model.matrix(~ 0 + meta_df$subtype)
  colnames(design) <- levels(meta_df$subtype)
  fit <- lmFit(expr_mat[members, , drop = FALSE], design)
  cont <- makeContrasts(contrasts = contrast_vec, levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  de <- topTable(fit2, number = Inf, sort.by = "P")
  de <- de %>% rownames_to_column("gene")
  return(head(de, top_n))
}

# 4C. Combined ranking: merge cor and de, compute mean rank
combine_leading <- function(cor_df, de_df){
  if(is.null(cor_df) && is.null(de_df)) return(NULL)
  if(is.null(cor_df)) cor_df <- data.frame(gene = de_df$gene, cor = NA)
  if(is.null(de_df)) de_df <- data.frame(gene = cor_df$gene, logFC = NA, t = NA, P.Value = NA)
  merged <- full_join(cor_df, de_df, by = "gene")
  merged <- merged %>%
    mutate(rank_cor = rank(-abs(cor), na.last = "keep"),
           rank_de  = rank(P.Value, na.last = "keep")) %>%
    mutate(mean_rank = rowMeans(cbind(ifelse(is.na(rank_cor), Inf, rank_cor),
                                      ifelse(is.na(rank_de), Inf, rank_de)), na.rm = TRUE)) %>%
    arrange(mean_rank)
  return(merged)
}

# ---------------------------
# 5. Run leading-edge extraction for each pathway in each subtype list
# ---------------------------
# prepare contrast strings for limma (we will compute all pairwise contrasts)
# ensure meta$subtype levels include TypeA, TypeB, TypeC
meta$subtype <- factor(meta$subtype, levels = c("TypeA","TypeB","TypeC"))
contrast_list <- list(
  "TypeA_vs_TypeB" = "TypeA - TypeB",
  "TypeA_vs_TypeC" = "TypeA - TypeC",
  "TypeB_vs_TypeC" = "TypeB - TypeC",
  "TypeB_vs_TypeA" = "TypeB - TypeA",
  "TypeC_vs_TypeA" = "TypeC - TypeA",
  "TypeC_vs_TypeB" = "TypeC - TypeB"
)

# helper to process a vector of pathways and return leading-edge tables
process_pathway_list <- function(pathways, subtype_label){
  res_list <- list()
  for(pw in pathways){
    cor_df <- leading_edge_cor(pw, gsva, expr, top_n = top_genes_per_pathway)
    # DE-based: use contrast where subtype is higher vs others
    # choose contrast depending on subtype_label
    if(subtype_label == "TypeA"){
      # A > B and A > C -> use TypeA - TypeB and TypeA - TypeC; combine DE ranks by min p
      de1 <- if(!is.null(gsets) && pw %in% names(gsets)) leading_edge_de(gsets[[pw]], expr, meta, "TypeA - TypeB", top_n = top_genes_per_pathway) else NULL
      de2 <- if(!is.null(gsets) && pw %in% names(gsets)) leading_edge_de(gsets[[pw]], expr, meta, "TypeA - TypeC", top_n = top_genes_per_pathway) else NULL
      # merge de1 and de2 by gene, take best p-value
      if(!is.null(de1) || !is.null(de2)){
        de_comb <- bind_rows(de1, de2) %>% group_by(gene) %>% summarize(logFC = mean(logFC, na.rm=TRUE), t = mean(t, na.rm=TRUE), P.Value = min(P.Value, na.rm=TRUE)) %>% arrange(P.Value)
      } else de_comb <- NULL
    } else if(subtype_label == "TypeB"){
      de1 <- if(!is.null(gsets) && pw %in% names(gsets)) leading_edge_de(gsets[[pw]], expr, meta, "TypeB - TypeA", top_n = top_genes_per_pathway) else NULL
      de2 <- if(!is.null(gsets) && pw %in% names(gsets)) leading_edge_de(gsets[[pw]], expr, meta, "TypeB - TypeC", top_n = top_genes_per_pathway) else NULL
      if(!is.null(de1) || !is.null(de2)){
        de_comb <- bind_rows(de1, de2) %>% group_by(gene) %>% summarize(logFC = mean(logFC, na.rm=TRUE), t = mean(t, na.rm=TRUE), P.Value = min(P.Value, na.rm=TRUE)) %>% arrange(P.Value)
      } else de_comb <- NULL
    } else { # TypeC
      de1 <- if(!is.null(gsets) && pw %in% names(gsets)) leading_edge_de(gsets[[pw]], expr, meta, "TypeC - TypeA", top_n = top_genes_per_pathway) else NULL
      de2 <- if(!is.null(gsets) && pw %in% names(gsets)) leading_edge_de(gsets[[pw]], expr, meta, "TypeC - TypeB", top_n = top_genes_per_pathway) else NULL
      if(!is.null(de1) || !is.null(de2)){
        de_comb <- bind_rows(de1, de2) %>% group_by(gene) %>% summarize(logFC = mean(logFC, na.rm=TRUE), t = mean(t, na.rm=TRUE), P.Value = min(P.Value, na.rm=TRUE)) %>% arrange(P.Value)
      } else de_comb <- NULL
    }
    
    combined <- combine_leading(cor_df, de_comb)
    # save per-pathway results
    res_list[[pw]] <- list(cor = cor_df, de = de_comb, combined = combined)
    # write CSVs
    if(!is.null(cor_df)) write.csv(cor_df, paste0(out_prefix, "_", subtype_label, "_", make.names(pw), "_cor_top.csv"), row.names = FALSE)
    if(!is.null(de_comb)) write.csv(de_comb, paste0(out_prefix, "_", subtype_label, "_", make.names(pw), "_de_top.csv"), row.names = FALSE)
    if(!is.null(combined)) write.csv(combined, paste0(out_prefix, "_", subtype_label, "_", make.names(pw), "_combined_top.csv"), row.names = FALSE)
  }
  return(res_list)
}

res_A <- process_pathway_list(pw_A, "TypeA")
res_B <- process_pathway_list(pw_B, "TypeB")
res_C <- process_pathway_list(pw_C, "TypeC")

# ---------------------------
# 6. Aggregate leading-edge genes into subtype signatures
# Strategy:
# - For each subtype, collect combined leading-edge gene lists across its pathways
# - Rank genes by: frequency of appearance across pathways, mean_rank (from combined), and mean absolute correlation
# - Produce final signature of top 'signature_size' genes
# ---------------------------
aggregate_signature <- function(res_list, expr_mat, gsva_mat){
  # collect combined tables
  combined_tables <- lapply(res_list, function(x) x$combined)
  combined_tables <- combined_tables[!sapply(combined_tables, is.null)]
  if(length(combined_tables) == 0) return(NULL)
  # stack with pathway label
  stacked <- bind_rows(lapply(names(combined_tables), function(pw){
    df <- combined_tables[[pw]]
    if(is.null(df)) return(NULL)
    df <- df %>% mutate(pathway = pw)
    return(df)
  }))
  if(is.null(stacked) || nrow(stacked) == 0) return(NULL)
  # frequency and mean_rank
  summary_df <- stacked %>%
    group_by(gene) %>%
    summarize(freq = n(),
              mean_rank = mean(mean_rank, na.rm = TRUE),
              mean_cor = mean(cor, na.rm = TRUE)) %>%
    arrange(desc(freq), mean_rank)
  # compute final score (lower is better): mean_rank / freq
  summary_df <- summary_df %>% mutate(final_score = mean_rank / freq) %>% arrange(final_score)
  # add expression-based stats
  summary_df$mean_expr <- rowMeans(expr_mat[intersect(summary_df$gene, rownames(expr_mat)), , drop = FALSE], na.rm = TRUE)
  return(summary_df)
}

sig_A_df <- aggregate_signature(res_A, expr, gsva)
sig_B_df <- aggregate_signature(res_B, expr, gsva)
sig_C_df <- aggregate_signature(res_C, expr, gsva)

# Save signatures (top N)
if(!is.null(sig_A_df)) write.csv(sig_A_df, paste0(out_prefix, "_TypeA_signature_allgenes.csv"), row.names = FALSE)
if(!is.null(sig_B_df)) write.csv(sig_B_df, paste0(out_prefix, "_TypeB_signature_allgenes.csv"), row.names = FALSE)
if(!is.null(sig_C_df)) write.csv(sig_C_df, paste0(out_prefix, "_TypeC_signature_allgenes.csv"), row.names = FALSE)

sigA_top <- if(!is.null(sig_A_df)) head(sig_A_df$gene, signature_size) else character(0)
sigB_top <- if(!is.null(sig_B_df)) head(sig_B_df$gene, signature_size) else character(0)
sigC_top <- if(!is.null(sig_C_df)) head(sig_C_df$gene, signature_size) else character(0)

write.csv(data.frame(gene = sigA_top), paste0(out_prefix, "_TypeA_signature_top", signature_size, ".csv"), row.names = FALSE)
write.csv(data.frame(gene = sigB_top), paste0(out_prefix, "_TypeB_signature_top", signature_size, ".csv"), row.names = FALSE)
write.csv(data.frame(gene = sigC_top), paste0(out_prefix, "_TypeC_signature_top", signature_size, ".csv"), row.names = FALSE)

cat("Signatures saved. Top counts: A =", length(sigA_top), "B =", length(sigB_top), "C =", length(sigC_top), "\n")

# ---------------------------
# 7. Quick diagnostics and summary outputs
# ---------------------------
# 7A. Overlap between signatures
overlap_AB <- intersect(sigA_top, sigB_top)
overlap_AC <- intersect(sigA_top, sigC_top)
overlap_BC <- intersect(sigB_top, sigC_top)
cat("Signature overlaps: A∩B =", length(overlap_AB), "A∩C =", length(overlap_AC), "B∩C =", length(overlap_BC), "\n")

# 7B. Save a summary report
report <- list(
  pathways = list(TypeA = pw_A, TypeB = pw_B, TypeC = pw_C),
  signature_topN = list(TypeA = sigA_top, TypeB = sigB_top, TypeC = sigC_top),
  overlaps = list(AB = overlap_AB, AC = overlap_AC, BC = overlap_BC)
)
saveRDS(report, paste0(out_prefix, "_report.rds"))

cat("All done. Outputs prefixed with:", out_prefix, "\n")
