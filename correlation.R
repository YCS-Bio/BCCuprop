# Required packages
install_if_missing <- function(pkgs){
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if(length(to_install)) install.packages(to_install, dependencies = TRUE)
}
install_if_missing(c("tidyverse","pheatmap","RColorBrewer","patchwork"))

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(patchwork)

# 1. Read CSV
df <- read.csv("senstivity_BC_Hub.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 2. Identify columns
#    - first 8 columns: drug sensitivity
#    - last 5 columns: gene expression
ncol_total <- ncol(df)
drug_cols <- colnames(df)[1:8]
gene_cols <- colnames(df)[(ncol_total-4):ncol_total]

# Optional: print to confirm
message("Drugs: ", paste(drug_cols, collapse = ", "))
message("Genes: ", paste(gene_cols, collapse = ", "))

# 3. Prepare long table for plotting and correlation
#    Keep sample identifier if present (rownames or a column). If no sample id, create one.
if (!("sample" %in% colnames(df))) {
  df$sample <- seq_len(nrow(df))
}

# long format for drugs and genes
drugs_long <- df %>%
  select(sample, all_of(drug_cols)) %>%
  pivot_longer(cols = -sample, names_to = "Drug", values_to = "Sensitivity")

genes_long <- df %>%
  select(sample, all_of(gene_cols)) %>%
  pivot_longer(cols = -sample, names_to = "Gene", values_to = "Expression")

# join to get all drug-gene pairs per sample
long_all <- left_join(drugs_long, genes_long, by = "sample")

# 4. Compute correlations for each Drug x Gene pair
corr_list <- long_all %>%
  group_by(Drug, Gene) %>%
  summarize(
    n = sum(!is.na(Sensitivity) & !is.na(Expression)),
    r_spearman = ifelse(n >= 3, cor(Sensitivity, Expression, method = "spearman", use = "pairwise.complete.obs"), NA_real_),
    p_spearman = ifelse(n >= 3, cor.test(Sensitivity, Expression, method = "spearman")$p.value, NA_real_),
    r_pearson = ifelse(n >= 3, cor(Sensitivity, Expression, method = "pearson", use = "pairwise.complete.obs"), NA_real_),
    p_pearson = ifelse(n >= 3, cor.test(Sensitivity, Expression, method = "pearson")$p.value, NA_real_),
    .groups = "drop"
  )

# Adjust p-values (BH) for multiple testing (use spearman p as default)
corr_list <- corr_list %>%
  mutate(
    p_spearman_adj = p.adjust(p_spearman, method = "BH"),
    p_pearson_adj = p.adjust(p_pearson, method = "BH")
  )

# Save results table
write.csv(corr_list, "correlations_table.csv", row.names = FALSE)
message("Saved correlations_table.csv")

# 5. Create heatmap of Spearman r (rows = Genes, cols = Drugs)
heat_mat <- corr_list %>%
  select(Drug, Gene, r_spearman) %>%
  pivot_wider(names_from = Drug, values_from = r_spearman) %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

# Create significance matrix for stars
sig_mat <- corr_list %>%
  mutate(star = case_when(
    is.na(p_spearman) ~ "",
    p_spearman < 0.001 ~ "***",
    p_spearman < 0.01  ~ "**",
    p_spearman < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  select(Drug, Gene, star) %>%
  pivot_wider(names_from = Drug, values_from = star) %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

# Heatmap colors
cols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(50)

# Plot heatmap and annotate with stars
pdf("cor_heatmap.pdf", width = 8, height = 6)
pheatmap(
  heat_mat,
  color = cols,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = sig_mat,       # show significance stars
  number_color = "black",
  fontsize_number = 14,
  main = "Spearman correlation (r) between Genes and Drugs",
  na_col = "grey90"
)
dev.off()
message("Saved cor_heatmap.pdf")

# 6. Scatter grid: one panel per Drug x Gene (40 plots)
#    We'll compute annotation text (r and p) and then facet by Gene ~ Drug
annot <- corr_list %>%
  select(Drug, Gene, r_spearman, p_spearman) %>%
  mutate(
    r_label = sprintf("r=%.2f", r_spearman),
    p_label = ifelse(is.na(p_spearman), "p=NA", sprintf("p=%.3g", p_spearman)),
    annot = paste(r_label, p_label, sep = "\n")
  )

plot_df <- left_join(long_all, annot, by = c("Drug", "Gene"))

# Basic scatter with linear fit and annotation; facet grid Gene rows x Drug columns
scatter_plot <- ggplot(plot_df, aes(x = Expression, y = Sensitivity)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.4) +
  facet_grid(rows = vars(Gene), cols = vars(Drug), scales = "free") +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  labs(x = "Expression", y = "Sensitivity") +
  geom_text(data = distinct(plot_df, Drug, Gene, annot), 
            aes(x = -Inf, y = Inf, label = annot),
            hjust = -0.1, vjust = 1.1, size = 3, inherit.aes = FALSE)

# Save scatter grid (large)
pdf("cor_scatter_grid.pdf", width = 16, height = 10)
print(scatter_plot)
dev.off()
message("Saved cor_scatter_grid.pdf")

# 7. Quick console summary
print(corr_list %>% arrange(p_spearman) %>% head(10))