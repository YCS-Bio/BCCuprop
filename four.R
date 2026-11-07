###############################################################################
# REO-based multiclass pipeline
# Models: svmLinear, ranger, glmnet, xgbTree
# Steps: training (nested CV), final training, test prediction, diagnostics,
#        survival analysis
###############################################################################

############################
# Libraries
############################
required_pkgs <- c(
  "data.table","dplyr","caret","pROC","xgboost",
  "ranger","glmnet","survival","survminer",
  "ggplot2","reshape2","scales","parallel","doParallel","RhpcBLASctl"
)
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)
}

library(data.table); library(dplyr); library(caret); library(pROC)
library(xgboost); library(ranger); library(glmnet)
library(survival); library(survminer)
library(ggplot2); library(reshape2); library(scales)
library(parallel); library(doParallel)

set.seed(2026)

############################
# User paths (EDIT)
############################
expr_train_file <- "expr_crg19_train.txt"   # genes x samples, first col gene
labels_train_file <- "labels_train.txt"     # sample, class
expr_test_file  <- "expr_crg19_test.txt"    # genes x samples, first col gene
time_file       <- "time.txt"               # sample, fustat, futime
outdir <- "model_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################
# Parameters
############################
top_k_genes <- 19
max_pairs   <- 100
outer_k     <- 5
nboot_stability <- 100
drop_frac   <- 0.10
flip_frac   <- 0.05
subsample_frac <- 0.8

############################
# Threading: all cores - 1
############################
detected_cores <- parallel::detectCores(logical = TRUE)
nthreads <- ifelse(is.na(detected_cores) || detected_cores <= 1, 1, detected_cores - 1)

Sys.setenv(OMP_NUM_THREADS = nthreads)
Sys.setenv(MKL_NUM_THREADS = nthreads)
Sys.setenv(OPENBLAS_NUM_THREADS = nthreads)
options(mc.cores = nthreads)
data.table::setDTthreads(nthreads)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(nthreads), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(nthreads), silent = TRUE)
}

cl <- makeCluster(nthreads)
registerDoParallel(cl)
message("Using ", nthreads, " threads (of ", detected_cores, " cores).")

############################
# Helper functions
############################
read_expr <- function(path){
  dt <- fread(path, data.table = FALSE)
  rn <- as.character(dt[[1]])
  mat <- as.data.frame(dt[,-1, drop = FALSE])
  rownames(mat) <- rn
  as.matrix(mat)
}

read_labels <- function(path){
  df <- fread(path, data.table = FALSE)
  colnames(df)[1:2] <- c("sample","class")
  df$class <- factor(df$class)
  df
}

select_top_pairs <- function(expr_mat, labels, top_k_genes = 200, max_pairs = 2000){
  # clean gene names
  genes <- rownames(expr_mat)
  genes <- trimws(genes)
  genes <- gsub("\\.\\d+$", "", genes)
  rownames(expr_mat) <- genes
  
  # remove duplicates
  if (any(duplicated(rownames(expr_mat)))) {
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  }
  
  gene_vars <- apply(expr_mat, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_k_genes, length(gene_vars))]
  
  combs <- t(combn(top_genes, 2))
  pairs <- paste0(combs[,1], "_gt_", combs[,2])
  scores <- rep(NA_real_, length(pairs))
  
  for (i in seq_along(pairs)){
    g1 <- combs[i,1]; g2 <- combs[i,2]
    
    if (!(g1 %in% rownames(expr_mat)) || !(g2 %in% rownames(expr_mat))) {
      scores[i] <- 1
      next
    }
    
    v <- as.integer(expr_mat[g1, ] > expr_mat[g2, ])
    tbl <- table(v, labels)
    pval <- tryCatch(chisq.test(tbl)$p.value, error = function(e) 1)
    scores[i] <- pval
  }
  
  ord <- order(scores, na.last = TRUE)
  pairs[ord][1:min(max_pairs, length(pairs))]
}

compute_reo_matrix <- function(expr_mat, pairs){
  pair_list <- strsplit(pairs, "_gt_")
  samples <- colnames(expr_mat)
  M <- matrix(NA_integer_, nrow = length(samples), ncol = length(pairs))
  colnames(M) <- pairs; rownames(M) <- samples
  
  for (j in seq_along(pair_list)){
    g1 <- pair_list[[j]][1]; g2 <- pair_list[[j]][2]
    if (g1 %in% rownames(expr_mat) && g2 %in% rownames(expr_mat)){
      M[,j] <- as.integer(expr_mat[g1, ] > expr_mat[g2, ])
    }
  }
  as.data.frame(M)
}

predict_probs_for_model <- function(model_obj, model_type, Xnew, class_levels){
  p <- predict(model_obj, Xnew, type = "prob")
  as.data.frame(p)[, class_levels, drop = FALSE]
}

eval_model_probs <- function(probs_mat, y_true){
  mat <- as.matrix(probs_mat)
  pred <- factor(colnames(mat)[max.col(mat, ties.method = "first")], levels = colnames(mat))
  acc <- mean(pred == y_true, na.rm = TRUE)
  
  per_class_auc <- sapply(colnames(mat), function(cl){
    y_bin <- as.integer(y_true == cl)
    if (length(unique(y_bin)) < 2) return(NA_real_)
    as.numeric(roc(y_bin, mat[,cl], quiet = TRUE)$auc)
  })
  mean_auc <- mean(per_class_auc, na.rm = TRUE)
  
  list(accuracy = acc, mean_ovr_auc = mean_auc)
}

############################
# Load data
############################
expr_train <- read_expr(expr_train_file)
expr_test  <- read_expr(expr_test_file)
labels_train <- read_labels(labels_train_file)

# align training samples
expr_train <- expr_train[, labels_train$sample, drop = FALSE]
y_train <- labels_train$class
class_levels <- levels(y_train)

############################
# Select REO pairs and build REO matrices
############################
selected_pairs <- select_top_pairs(expr_train, y_train, top_k_genes, max_pairs)
writeLines(selected_pairs, file.path(outdir, "selected_pairs.txt"))

X_train <- compute_reo_matrix(expr_train, selected_pairs)
X_test  <- compute_reo_matrix(expr_test,  selected_pairs)

# drop columns that are all NA in train or test
keep_idx <- which(!(apply(X_train, 2, function(x) all(is.na(x))) |
                      apply(X_test,  2, function(x) all(is.na(x)))))
X_train <- X_train[, keep_idx, drop = FALSE]
X_test  <- X_test[,  keep_idx, drop = FALSE]
selected_pairs <- selected_pairs[keep_idx]

############################
# Nested CV (svm, xgb, rf, enet)
############################
outer_folds <- createFolds(y_train, k = outer_k)
model_types <- c("svm","xgb","rf","enet")
outer_results <- vector("list", outer_k)

ctrl_inner <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  allowParallel = TRUE
)

for (i in seq_along(outer_folds)){
  val_idx <- outer_folds[[i]]
  tr_idx  <- setdiff(seq_len(nrow(X_train)), val_idx)
  
  Xtr <- X_train[tr_idx, , drop = FALSE]
  Xval <- X_train[val_idx, , drop = FALSE]
  ytr <- y_train[tr_idx]
  yval <- y_train[val_idx]
  
  eval_list <- list()
  
  # SVM
  svm_model <- tryCatch(
    train(Xtr, ytr, method = "svmLinear", trControl = ctrl_inner,
          tuneGrid = expand.grid(C = c(0.01, 0.1, 1, 10))),
    error = function(e) NULL
  )
  if (!is.null(svm_model))
    eval_list$svm <- eval_model_probs(predict(svm_model, Xval, type = "prob"), yval)
  
  # XGB
  xgb_model <- tryCatch(
    train(Xtr, ytr, method = "xgbTree", trControl = ctrl_inner, tuneLength = 5),
    error = function(e) NULL
  )
  if (!is.null(xgb_model))
    eval_list$xgb <- eval_model_probs(predict(xgb_model, Xval, type = "prob"), yval)
  
  # RF
  rf_model <- tryCatch(
    train(Xtr, ytr, method = "ranger", trControl = ctrl_inner,
          tuneLength = 3, num.threads = nthreads),
    error = function(e) NULL
  )
  if (!is.null(rf_model))
    eval_list$rf <- eval_model_probs(predict(rf_model, Xval, type = "prob"), yval)
  
  # ENet
  enet_model <- tryCatch(
    train(Xtr, ytr, method = "glmnet", trControl = ctrl_inner, tuneLength = 5),
    error = function(e) NULL
  )
  if (!is.null(enet_model))
    eval_list$enet <- eval_model_probs(predict(enet_model, Xval, type = "prob"), yval)
  
  outer_results[[i]] <- eval_list
  cat("Outer fold", i, "done.\n")
}

############################
# Select best model family (accuracy)
############################
mean_acc <- sapply(model_types, function(mt){
  vals <- sapply(outer_results, function(o) if (!is.null(o[[mt]])) o[[mt]]$accuracy else NA_real_)
  mean(vals, na.rm = TRUE)
})
print(mean_acc)
best_model_type <- names(which.max(mean_acc))
cat("Best model type:", best_model_type, "\n")

############################
# Final training (svm, xgb, rf, enet)
############################
final_models <- list()

ctrl_final <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  allowParallel = TRUE
)

final_models$svm <- tryCatch(
  train(X_train, y_train, method = "svmLinear", trControl = ctrl_final,
        tuneGrid = expand.grid(C = c(0.01, 0.1, 1, 10))),
  error = function(e) NULL
)

final_models$xgb <- tryCatch(
  train(X_train, y_train, method = "xgbTree", trControl = ctrl_final, tuneLength = 5),
  error = function(e) NULL
)

final_models$rf <- tryCatch(
  train(X_train, y_train, method = "ranger", trControl = ctrl_final,
        tuneLength = 3, num.threads = nthreads),
  error = function(e) NULL
)

final_models$enet <- tryCatch(
  train(X_train, y_train, method = "glmnet", trControl = ctrl_final, tuneLength = 5),
  error = function(e) NULL
)

saveRDS(final_models, file.path(outdir, "final_models.rds"))

############################
# Predict on test (no test labels)
############################
probs_list <- list()
for (mt in names(final_models)){
  if (!is.null(final_models[[mt]])){
    probs_list[[mt]] <- predict_probs_for_model(final_models[[mt]], mt, X_test, class_levels)
  }
}
if (length(probs_list) == 0) stop("No final models available.")

ensemble_probs_test <- Reduce("+", probs_list) / length(probs_list)
rownames(ensemble_probs_test) <- rownames(X_test)
predicted_class <- colnames(ensemble_probs_test)[max.col(ensemble_probs_test, ties.method = "first")]

write.table(
  data.frame(sample = rownames(ensemble_probs_test), predicted = predicted_class),
  file.path(outdir, "test_predictions.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

############################
# Diagnostics: distribution, stability, ensemble agreement
############################
# Train-set ensemble probabilities (optimistic reference)
probs_train_list <- list()
for (mt in names(final_models)){
  if (!is.null(final_models[[mt]])){
    probs_train_list[[mt]] <- predict_probs_for_model(final_models[[mt]], mt, X_train, class_levels)
  }
}
ensemble_probs_train <- Reduce("+", probs_train_list) / length(probs_train_list)
rownames(ensemble_probs_train) <- rownames(X_train)

# Distribution comparison (train vs test)
summary_per_class <- data.frame(
  class = class_levels,
  train_mean = sapply(class_levels, function(cl) mean(ensemble_probs_train[,cl], na.rm = TRUE)),
  test_mean  = sapply(class_levels, function(cl) mean(ensemble_probs_test[,cl],  na.rm = TRUE))
)
write.table(summary_per_class, file.path(outdir, "prob_dist_train_vs_test_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

ks_res <- lapply(class_levels, function(cl){
  tr <- ensemble_probs_train[,cl]; te <- ensemble_probs_test[,cl]
  ks <- tryCatch(ks.test(tr, te), error = function(e) NULL)
  data.frame(
    class = cl,
    ks_stat = ifelse(is.null(ks), NA, ks$statistic),
    ks_pvalue = ifelse(is.null(ks), NA, ks$p.value)
  )
})
ks_table <- do.call(rbind, ks_res)
write.table(ks_table, file.path(outdir, "ks_train_vs_test_per_class.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

m_tr <- melt(as.data.frame(ensemble_probs_train), variable.name = "class", value.name = "prob")
m_tr$set <- "train"
m_te <- melt(as.data.frame(ensemble_probs_test),  variable.name = "class", value.name = "prob")
m_te$set <- "test"
m_all <- rbind(m_tr, m_te)
p_density <- ggplot(m_all, aes(x = prob, color = set, fill = set)) +
  geom_density(alpha = 0.15) +
  facet_wrap(~ class, ncol = 2) +
  theme_minimal()
ggsave(file.path(outdir, "prob_density_train_vs_test.png"), p_density, width = 10, height = 6, dpi = 300)

# Stability: perturb test REO matrix
Xmat <- as.matrix(X_test)
n_samples <- nrow(Xmat)
orig_pred <- predicted_class

ensemble_predict_from_X <- function(Xmat_local){
  probs_list_b <- list()
  for (mt in names(final_models)){
    if (!is.null(final_models[[mt]])){
      probs_list_b[[mt]] <- predict_probs_for_model(final_models[[mt]], mt, as.data.frame(Xmat_local), class_levels)
    }
  }
  ens <- Reduce("+", probs_list_b) / length(probs_list_b)
  colnames(ens)[max.col(ens, ties.method = "first")]
}

# A) dropout
preds_drop <- matrix(NA_character_, nrow = n_samples, ncol = nboot_stability)
for (b in 1:nboot_stability){
  drop_idx <- sample(seq_len(ncol(Xmat)), size = max(1, floor(drop_frac * ncol(Xmat))))
  Xpert <- Xmat
  Xpert[, drop_idx] <- Xpert[sample(nrow(Xpert)), drop_idx]
  preds_drop[,b] <- ensemble_predict_from_X(Xpert)
}
stab_drop <- rowMeans(preds_drop == orig_pred)

# B) bit flips
preds_flip <- matrix(NA_character_, nrow = n_samples, ncol = nboot_stability)
for (b in 1:nboot_stability){
  Xpert <- Xmat
  nflip <- max(1, floor(flip_frac * ncol(Xmat) * nrow(Xmat)))
  rows <- sample(seq_len(nrow(Xpert)), nflip, replace = TRUE)
  cols <- sample(seq_len(ncol(Xpert)), nflip, replace = TRUE)
  for (k in seq_along(rows)){
    Xpert[rows[k], cols[k]] <- 1L - Xpert[rows[k], cols[k]]
  }
  preds_flip[,b] <- ensemble_predict_from_X(Xpert)
}
stab_flip <- rowMeans(preds_flip == orig_pred)

# C) subsample features
# C) subsample features (XGBoost-safe version)
preds_sub <- matrix(NA_character_, nrow = n_samples, ncol = nboot_stability)

for (b in 1:nboot_stability) {
  
  # choose subset of features to KEEP unperturbed
  keep_idx <- sample(seq_len(ncol(Xmat)),
                     size = max(2, floor(subsample_frac * ncol(Xmat))))
  
  # start with full matrix
  Xpert <- Xmat
  
  # perturb ONLY the features NOT in keep_idx
  drop_idx <- setdiff(seq_len(ncol(Xmat)), keep_idx)
  
  if (length(drop_idx) > 0) {
    # randomize these columns (preserves 0/1 structure)
    Xpert[, drop_idx] <- Xpert[sample(nrow(Xpert)), drop_idx]
  }
  
  # predict using full feature set (safe for XGBoost)
  preds_sub[, b] <- ensemble_predict_from_X(Xpert)
}

stab_sub <- rowMeans(preds_sub == orig_pred, na.rm = TRUE)

stability_df <- data.frame(
  sample = rownames(X_test),
  predicted = orig_pred,
  stab_drop = stab_drop,
  stab_flip = stab_flip,
  stab_sub  = stab_sub
)
write.table(stability_df, file.path(outdir, "prediction_stability.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Ensemble agreement
preds_per_model <- as.data.frame(lapply(final_models, function(m){
  if (is.null(m)) return(NULL)
  p <- predict_probs_for_model(m, "dummy", X_test, class_levels)
  colnames(p)[max.col(p, ties.method = "first")]
}))
preds_per_model <- preds_per_model[, !sapply(preds_per_model, is.null), drop = FALSE]
rownames(preds_per_model) <- rownames(X_test)

majority_vote <- apply(preds_per_model, 1, function(x){
  tab <- table(x); names(tab)[which.max(tab)]
})
agreement_fraction <- sapply(seq_len(nrow(preds_per_model)), function(i){
  mean(preds_per_model[i, ] == majority[i])
})

ensemble_agree_df <- data.frame(
  sample = rownames(preds_per_model),
  majority = majority_vote,
  agreement = agreement_fraction
)
write.table(ensemble_agree_df, file.path(outdir, "ensemble_agreement.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Confidence metrics
conf_df <- data.frame(
  sample = rownames(ensemble_probs_test),
  predicted = orig_pred,
  max_prob = apply(ensemble_probs_test, 1, max),
  margin   = apply(ensemble_probs_test, 1, function(x){
    s <- sort(x, decreasing = TRUE); s[1] - s[2]
  }),
  stability_drop = stab_drop,
  stability_flip = stab_flip,
  stability_sub  = stab_sub,
  agreement      = agreement_fraction
)
write.table(conf_df, file.path(outdir, "prediction_confidence.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

############################
# Survival analysis
############################
time_df <- fread(time_file, data.table = FALSE)
colnames(time_df)[1:3] <- c("sample","fustat","futime")

# Train KM (true labels)
train_meta <- merge(
  data.frame(sample = labels_train$sample, cluster = labels_train$class),
  time_df, by = "sample", all.x = TRUE
)
if (all(c("futime","fustat") %in% colnames(train_meta))){
  fit_train <- survfit(Surv(futime, fustat) ~ cluster, data = train_meta)
  p_train <- ggsurvplot(fit_train, data = train_meta, pval = TRUE)
  ggsave(file.path(outdir, "KM_train_clusters.png"), p_train$plot,
         width = 6, height = 5, dpi = 300)
}

# Test KM (predicted clusters)
test_meta <- merge(
  data.frame(sample = rownames(ensemble_probs_test), predicted = predicted_class),
  time_df, by = "sample", all.x = TRUE
)
if (all(c("futime","fustat") %in% colnames(test_meta))){
  fit_test <- survfit(Surv(futime, fustat) ~ predicted, data = test_meta)
  p_test <- ggsurvplot(fit_test, data = test_meta, pval = TRUE)
  ggsave(file.path(outdir, "KM_test_predicted_clusters.png"), p_test$plot,
         width = 6, height = 5, dpi = 300)
}

############################
# Save artifacts and cleanup
############################
saveRDS(list(
  final_models = final_models,
  outer_results = outer_results,
  selected_pairs = selected_pairs,
  ensemble_probs_test = ensemble_probs_test,
  conf_df = conf_df,
  stability = stability_df,
  agreement = ensemble_agree_df
), file.path(outdir, "pipeline_results.rds"))

stopCluster(cl)
registerDoSEQ()
message("Pipeline complete. Outputs written to: ", outdir)
