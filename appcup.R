library(shiny)
library(tidyverse)
library(pheatmap)
library(ggplot2)

# ---- Load models and preprocessing ----
final_enet <- readRDS("final_enet.rds")
final_rf   <- readRDS("final_rf.rds")
preproc    <- readRDS("preprocessing_artifacts.rds")

# ---- Helper: preprocess expression matrix ----
preprocess_expr <- function(expr_mat, preproc) {
  common_genes <- intersect(preproc$genes, rownames(expr_mat))
  expr_sub <- expr_mat[common_genes, , drop = FALSE]
  
  means_sub <- preproc$gene_means[common_genes]
  sds_sub   <- preproc$gene_sds[common_genes]
  
  expr_scaled <- sweep(expr_sub, 1, means_sub, "-")
  expr_scaled <- sweep(expr_scaled, 1, sds_sub, "/")
  
  t(expr_scaled)  # samples x genes
}

# ---- Demo dataset ----
demo_expr <- function() {
  set.seed(1)
  mat <- matrix(rnorm(19 * 5, 6, 1), nrow = 19, ncol = 5)
  rownames(mat) <- preproc$genes
  colnames(mat) <- paste0("Demo_", 1:5)
  mat
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Cuproptosis 19‑Gene Classifier"),
  
  tabsetPanel(
    
    # ---------------- TAB 1: Upload & QC ----------------
    tabPanel("Upload & QC",
             sidebarLayout(
               sidebarPanel(
                 fileInput("expr_file", "Upload expression matrix (CSV)",
                           accept = c(".csv")),
                 actionButton("demo_btn", "Load Demo Dataset"),
                 br(), br(),
                 helpText(
                   strong("File requirements:"),
                   br(),
                   "• CSV file with genes as rows and samples as columns.",
                   br(),
                   "• Gene names must be HGNC symbols.",
                   br(),
                   "• Expression values must be log2-transformed.",
                   br(),
                   "• Supports one or multiple patients.",
                   br(),
                   "• Missing genes will be shown below."
                 )
               ),
               mainPanel(
                 h4("QC: Missing Genes"),
                 tableOutput("qc_table")
               )
             )
    ),
    
    # ---------------- TAB 2: Predictions ----------------
    tabPanel("Predictions",
             sidebarLayout(
               sidebarPanel(
                 actionButton("run_btn", "Run Prediction"),
                 br(), br(),
                 downloadButton("download_pred", "Download Predictions")
               ),
               mainPanel(
                 h4("Prediction Results"),
                 tableOutput("pred_table")
               )
             )
    ),
    
    # ---------------- TAB 3: Heatmap ----------------
    tabPanel("Heatmap",
             h4("19‑Gene Expression Heatmap"),
             plotOutput("heatmap_plot", height = "600px")
    ),
    
    # ---------------- TAB 4: PCA ----------------
    tabPanel("PCA",
             h4("PCA of Uploaded Samples"),
             plotOutput("pca_plot", height = "600px")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  expr_data <- reactiveVal(NULL)
  predictions <- reactiveVal(NULL)
  
  # ---- Load demo dataset ----
  observeEvent(input$demo_btn, {
    expr_data(demo_expr())
    output$qc_table <- renderTable({
      tibble(Message = "All 19 model genes are present — no missing genes.")
    })
  })
  
  # ---- Load uploaded file ----
  observeEvent(input$expr_file, {
    req(input$expr_file)
    expr <- read.csv(input$expr_file$datapath, row.names = 1, check.names = FALSE)
    expr_data(as.matrix(expr))
    
    # QC
    missing_genes <- setdiff(preproc$genes, rownames(expr))
    output$qc_table <- renderTable({
      if (length(missing_genes) == 0) {
        tibble(Message = "All 19 model genes are present — no missing genes.")
      } else {
        tibble(Missing_Genes = missing_genes)
      }
    })
  })
  
  # ---- Run prediction ----
  observeEvent(input$run_btn, {
    req(expr_data())
    expr <- expr_data()
    
    # Preprocess
    X <- preprocess_expr(expr, preproc)
    
    # Predict
    prob_enet <- predict(final_enet, X, type = "prob")[, 2]
    prob_rf   <- predict(final_rf,   X, type = "prob")[, 2]
    
    class_enet <- ifelse(prob_enet >= 0.5, "Cluster 2", "Cluster 1")
    class_rf   <- ifelse(prob_rf   >= 0.5, "Cluster 2", "Cluster 1")
    
    consensus <- ifelse(
      (class_enet == "Cluster 2") + (class_rf == "Cluster 2") >= 1,
      "Cluster 2", "Cluster 1"
    )
    
    conf_enet <- abs(prob_enet - 0.5) * 2
    conf_rf   <- abs(prob_rf   - 0.5) * 2
    conf_consensus <- (conf_enet + conf_rf) / 2
    
    pred_df <- tibble(
      Sample          = rownames(X),
      ENET_Prob       = round(prob_enet, 4),
      ENET_Class      = class_enet,
      RF_Prob         = round(prob_rf, 4),
      RF_Class        = class_rf,
      Consensus_Class = consensus,
      Confidence      = round(conf_consensus, 4)
    )
    
    predictions(pred_df)
    
    output$pred_table <- renderTable(pred_df)
  })
  
  # ---- Download predictions ----
  output$download_pred <- downloadHandler(
    filename = function() {
      paste0("cuproptosis_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(predictions(), file, row.names = FALSE)
    }
  )
  
  # ---- Heatmap ----
  output$heatmap_plot <- renderPlot({
    req(expr_data())
    expr <- expr_data()
    common <- intersect(preproc$genes, rownames(expr))
    mat <- expr[common, , drop = FALSE]
    pheatmap(mat, scale = "row", clustering_distance_rows = "correlation")
  })
  
  # ---- PCA ----
  output$pca_plot <- renderPlot({
    req(expr_data(), predictions())
    expr <- expr_data()
    common <- intersect(preproc$genes, rownames(expr))
    mat <- t(expr[common, , drop = FALSE])
    
    pca <- prcomp(mat, scale. = TRUE)
    df <- data.frame(
      PC1 = pca$x[,1],
      PC2 = pca$x[,2],
      Sample = rownames(mat),
      Cluster = predictions()$Consensus_Class
    )
    
    ggplot(df, aes(PC1, PC2, color = Cluster, label = Sample)) +
      geom_point(size = 4) +
      geom_text(vjust = -1) +
      theme_minimal() +
      scale_color_manual(values = c("Cluster 1" = "#1b9e77", "Cluster 2" = "#d95f02"))
  })
}

shinyApp(ui, server)