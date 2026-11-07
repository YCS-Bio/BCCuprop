library(tidyverse)
library(forcats)
library(patchwork)

# Load data
clinical <- read.csv("clinical_filtered_breast_TNM.csv", stringsAsFactors = FALSE)

# Filter tumor samples
clinical <- clinical %>% 
  filter(tissue_type == "Tumor")

# Create age category
clinical <- clinical %>%
  mutate(
    age_group = case_when(
      age_dx < 40 ~ "<40",
      age_dx >= 40 & age_dx < 50 ~ "40-50",
      age_dx >= 50 & age_dx < 60 ~ "50-60",
      age_dx >= 60 & age_dx < 70 ~ "60-70",
      age_dx >= 70 ~ ">70",
      TRUE ~ NA_character_
    )
  )

# Variables to plot
vars <- c("gender", "race", "ethnicity", "age_group", "stage")

# Minimal palette reused across variables
min_colors <- c("#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD")

# Function to create one stacked bar per variable
make_bar <- function(df, var) {
  
  df <- df %>% mutate(category = .data[[var]])
  
  # Apply custom ordering rules
  if (var == "race") {
    df$category <- fct_explicit_na(df$category, na_level = "Not Reported")
    df$category <- fct_relevel(
      df$category,
      "American Indian or Alaska native",
      "asian",
      "black or african american",
      "white",
      "Not Reported",
      after = Inf
    )
  }
  
  if (var == "ethnicity") {
    df$category <- fct_explicit_na(df$category, na_level = "Not Reported")
    df$category <- fct_collapse(
      df$category,
      "Not Reported" = c("not reported", "unknown", "Not Reported")
    )
    df$category <- fct_relevel(
      df$category,
      "not hispanic or latino",
      "hispanic or latino",
      "Not Reported",
      after = Inf
    )
  }
  
  if (var == "age_group") {
    df$category <- fct_explicit_na(df$category, na_level = "Not Reported")
    df$category <- fct_relevel(
      df$category,
      "<40", "40-50", "50-60", ">60", "Not Reported",
      after = Inf
    )
  }
  
  # For other variables, keep NA at bottom
  if (!(var %in% c("race", "ethnicity", "age_group"))) {
    df$category <- fct_explicit_na(df$category, na_level = "Not Reported")
    df$category <- fct_relevel(df$category, "Not Reported", after = Inf)
  }
  
  df %>%
    count(category) %>%
    ggplot(aes(x = var, y = n, fill = category)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = n),
              position = position_stack(vjust = 0.5),
              size = 4) +
    scale_fill_manual(values = min_colors) +
    labs(x = "", y = "Count", fill = "", title = var) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# Build each bar
plots <- lapply(vars, function(v) make_bar(clinical, v))

# Combine into one row
combined_plot <- wrap_plots(plots, nrow = 1)

combined_plot
ggsave(
  filename = "clinical_stacked_bars.png",
  plot = combined_plot,
  width = 12,
  height = 5,
  dpi = 300
)
