# ------------------------------------------------------------
# Input:
#   expr:   gene expression matrix (genes x samples)
#   y:      binary vector (0/1) for anaerobic vs aerobic class
#   pairs:  data.frame with columns geneA, geneB
# Output:
#   ranked table with Fisher p, AUC, direction, and counts
# ------------------------------------------------------------

library(pROC)
library(dplyr)

screen_reo_pairs <- function(expr, y, pairs) {
  
  # Ensure gene names match
  rownames(expr) <- toupper(rownames(expr))
  pairs$geneA <- toupper(pairs$geneA)
  pairs$geneB <- toupper(pairs$geneB)
  
  results <- lapply(1:nrow(pairs), function(i) {
    gA <- pairs$geneA[i]
    gB <- pairs$geneB[i]
    
    if (!(gA %in% rownames(expr)) | !(gB %in% rownames(expr))) {
      return(data.frame(
        geneA = gA, geneB = gB,
        fisher_p = NA, auc = NA,
        direction = NA,
        n_A_gt_B = NA, n_B_gt_A = NA
      ))
    }
    
    # REO: A > B?
    reo <- expr[gA, ] > expr[gB, ]
    
    # --- Fisher's exact test (only if both reo and y have variation) ---
    if (length(unique(reo)) < 2 || length(unique(y)) < 2) {
      fisher_p <- NA
    } else {
      tab <- table(reo, y)
      if (nrow(tab) == 2 && ncol(tab) == 2) {
        fisher_p <- fisher.test(tab)$p.value
      } else {
        fisher_p <- NA
      }
    }
    
    # --- AUC: treat REO TRUE/FALSE as predictor ---
    auc_val <- tryCatch({
      if (length(unique(reo)) == 2) {
        as.numeric(pROC::roc(y, as.numeric(reo))$auc)
      } else {
        NA
      }
    }, error = function(e) NA)
    
    # --- Direction (requires both groups exist) ---
    direction <- NA
    if (length(unique(y)) == 2 && length(unique(reo)) == 2) {
      mean_aerobic <- mean(reo[y == 1], na.rm = TRUE)   # assuming 1 = aerobic
      mean_anaerobic <- mean(reo[y == 0], na.rm = TRUE) # assuming 0 = anaerobic
      if (!is.na(mean_aerobic) && !is.na(mean_anaerobic)) {
        direction <- ifelse(mean_aerobic > mean_anaerobic,
                            "A>B in aerobic",
                            "A>B in anaerobic")
      }
    }
    
    data.frame(
      geneA = gA,
      geneB = gB,
      fisher_p = fisher_p,
      auc = auc_val,
      direction = direction,
      n_A_gt_B = sum(reo, na.rm = TRUE),
      n_B_gt_A = sum(!reo, na.rm = TRUE)
    )
  })
  
  results <- do.call(rbind, results)
  
  # Rank by Fisher p then AUC
  results %>%
    dplyr::arrange(fisher_p, dplyr::desc(auc))
}
setwd("D:/aerobic")
expr=read.table("expr_aerobic.txt",row.names = 1, check.names = F, sep = "\t",header=T )
y=read.table("y.txt",row.names = 1, check.names = F, sep = "\t",header=T )
expr=expr[rownames(y)]
y <- y[colnames(expr), 1]   # 按 expr 列名排序，取第一列
y <- factor(y)                  # 转为因子（适合分类分析）
pairs <- data.frame(
  geneA = c(
    "NDUFS1","NDUFV1","NDUFA9","COX4I1","COX5A","UQCRC1","CYC1","ATP5F1A","ATP5F1B",
    "CS","IDH3A","OGDH","DLST","MDH2","SDHB",
    "COX4I1","COX4I1","COX4I1","ATP5F1A","ATP5F1A","CS","CS","NDUFS1","NDUFS1",
    "PDHA1","PDHB","DLAT","DLD","CS","ATP5F1A"
  ),
  geneB = c(
    "LDHA","LDHA","LDHA","LDHA","LDHA","LDHA","LDHA","LDHA","LDHA",
    "LDHA","LDHA","LDHA","LDHA","LDHA","LDHA",
    "HK2","PFKP","PKM","PFKP","PKM","PFKP","PKM","PFKP","PKM",
    "PDK1","PDK1","PDK1","PDK1","PDK1","PDK1"
  )
)
results <- screen_reo_pairs(expr, y, pairs)
head(results, 10)
