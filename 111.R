####====================data preparation====================####
rm(list=ls())
options(stringaAsFactors=F)
setwd("D:\\013GEO") 

library(caret)  # for stratified split

df_exp <- read.table("merged_normalized_data.txt", header=TRUE, sep="\t", check.names=FALSE)

# Quick structure check
str(df_exp)
summary(df_exp$futime)
table(df_exp$fustat)  # should be 0/1


# Ensure proper types
df_exp$futime  <- as.numeric(df_exp$futime)
df_exp$fustat  <- as.numeric(df_exp$fustat)
#df_exp$fustat  <- ifelse(df_exp$fustat > 1, 1, df_exp$fustat)  # force binary if needed

df_exp$age              <- as.numeric(df_exp$age)
df_exp$lymphnode      <- factor(df_exp$lymphnode)        # e.g., N0/N1/N2
df_exp$tstage <- factor(df_exp$tstage)   # e.g., T1/T2/T3/T4
df_exp$subtype    <- factor(df_exp$subtype)

# Remove impossible/failing rows
df_exp <- df_exp[is.finite(df_exp$futime) & df_exp$futime > 0, ]



####====================Train/validation splits====================####

library(caret)  # for stratified split
set.seed(5201314)


# Stratified split on event
idx <- createDataPartition(df_exp$fustat, p = 0.7, list = FALSE)
train <- df_exp[idx, ]
vali  <- df_exp[-idx, ]

list_train_vali_Data <- list(Dataset1 = train, Dataset2 = vali)

####====================Univariate Cox screening and basic baselines====================####
##- Purpose: Confirm signals and reduce noise before ML modeling.
##- Run univariate Cox for each gene and key covariates.
library(survival)
genelist <- c("ALB", "APC", "AQP1", "ARF1", "ATP7A", "ATP7B", "BACE1", "CCS", "CDK1", "CDKN2A", "CYP1A1", "DAXX", "DLD", "DLST", "LIPT1", "MTF1", "NFAT5", "NFE2L2", "PAM", "PARK7", "PDE3B", "PPP1R15A", "PRNP", "SCO1", "SCO2", "SNCA", "SNCG")

uni_cox <- function(dat, vars){
  out <- lapply(vars, function(v){
    f <- as.formula(paste0("Surv(futime, fustat) ~ `", v, "`"))
    m <- coxph(f, data = dat)
    s <- summary(m)
    c(var = v,
      HR = s$coefficients[1,"exp(coef)"],
      z  = s$coefficients[1,"z"],
      p  = s$coefficients[1,"Pr(>|z|)"])
  })
  as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
}

vars_to_test <- c(genelist, "age","lymphnode","tstage","subtype")
uni_res <- uni_cox(train, vars_to_test)
uni_res$p <- as.numeric(uni_res$p)
uni_sig <- subset(uni_res, p < 0.05)
uni_sig

## Interpretation: Keep genes/covariates with strong association; optionally use FDR correction.

uni_res$padj <- p.adjust(uni_res$p, method = "BH")
uni_sig_fdr <- subset(uni_res, padj < 0.05)

####====================Run Mime1’s 101-model ensemble====================####
## Prepare: Mime1 expects consistent dataframes and a gene list. You’ve already got log2 expression, factors set.
## Set options: Use univariate Cox filter inside Mime1; include clinical covariates if supported as features.
# 1) Define your candidate genes (27 provided)
genelist <- c("ALB","APC","AQP1","ARF1","ATP7A","ATP7B","BACE1","CCS","CDK1","CDKN2A",
              "CYP1A1","DAXX","DLD","DLST","LIPT1","MTF1","NFAT5","NFE2L2","PAM","PARK7",
              "PDE3B","PPP1R15A","PRNP","SCO1","SCO2","SNCA","SNCG")

# 2) Rename and coerce survival columns; create ID if missing
normalize_one <- function(d){
  d <- as.data.frame(d)
  
  # Create ID if not present
  if(!"ID" %in% colnames(d)){
    if (!is.null(rownames(d)) && all(rownames(d) != "")) {
      d$ID <- rownames(d)
    } else if ("Patient" %in% colnames(d)) {
      d$ID <- d$Patient
    } else {
      d$ID <- seq_len(nrow(d))
    }
  }
  
  # Map futime/fustat to required names if needed
  if("OS.time" %in% colnames(d) == FALSE && "futime" %in% colnames(d)) {
    names(d)[names(d) == "futime"] <- "OS.time"
  }
  if("OS" %in% colnames(d) == FALSE && "fustat" %in% colnames(d)) {
    names(d)[names(d) == "fustat"] <- "OS"
  }
  
  # Type coercions
  d$OS.time <- as.numeric(d$OS.time)
  d$OS      <- as.numeric(d$OS)
  d$OS[is.na(d$OS)] <- 0
  d$OS[d$OS > 1] <- 1  # enforce binary
  
  # Reorder columns: ID, OS.time, OS, then the rest
  first3 <- c("ID","OS.time","OS")
  rest   <- setdiff(colnames(d), first3)
  d <- d[, c(first3, rest)]
  
  # Remove rows with invalid time
  d <- d[is.finite(d$OS.time) & d$OS.time > 0, ]
  
  return(d)
}

list_train_vali_Data <- lapply(list_train_vali_Data, normalize_one)

# 3) Find genes present in ALL cohorts and restrict to them
common_genes <- Reduce(intersect, lapply(list_train_vali_Data, function(d) {
  intersect(genelist, colnames(d))
}))

# Warn if you lose genes
message(sprintf("Candidate genes requested: %d; common across cohorts: %d",
                length(genelist), length(common_genes)))

# 4) Drop non-common genes from each dataset to align features
list_train_vali_Data <- lapply(list_train_vali_Data, function(d) {
  keep_cols <- c("ID","OS.time","OS", common_genes)
  d[, keep_cols, drop = FALSE]
})
# z-score normalized already done

# 6) Choose training dataset
train_data <- list_train_vali_Data[[1]]  # e.g., Dataset1 after lapply
var.select <- randomForestSRC::var.select.rfsrc


library(GSEABase)
library(GSVA)
library(cancerclass)
library(mixOmics)
library(sparrow)
library(sva)
library(CoxBoost)
library(fastAdaboost)
library(randomForestSRC)
library(varPro)
library(Mime1)
out_file <- "my_ML.Dev.Prog.Sig_modified_functions.R"
source(out_file)
res <- ML.Dev.Prog.Sig_modified(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  unicox.filter.for.candi = TRUE,
  unicox_p_cutoff = 0.05,
  candidate_genes = common_genes,  # use harmonized gene set
  mode = "all", # exclude RSF
  nodesize = 5,
  seed = 5201314
)

# Check what components are inside res
str(res, max.level = 1)
names(res)

# Archive the entire result object for reproducibility
saveRDS(res, "res_full.rds")
#Export each component separatel
#1. Concordance index results (Cindex.res)
## a dataframe already
write.csv(res$Cindex.res, "Cindex_results.csv", row.names = FALSE)
#2. Machine learning results (ml.res)
##This is a list of 117 elements. You’ll likely want to save the whole list as RDS:
saveRDS(res$ml.res, "ml_results.rds")
#3. Risk scores (riskscore)
##Same structure as ml.res (list of 117). Save as RDS:
saveRDS(res$riskscore, "riskscores.rds")
#4. Signature genes (Sig.genes)
##This is a character vector of 14 genes:
write.csv(data.frame(Sig.genes = res$Sig.genes), "signature_genes.csv", row.names = FALSE)
  
#Save a summary log
sink("res_summary.txt")
print(summary(res))
sink()

####Save as an R object (for later reuse####
# Or save multiple objects together
save(res, train_data, list_train_vali_Data, file = "analysis_results.RData")

# Load back later
##load("res.RData")
# Reload all objects saved in the workspace
load("analysis_results.RData")




####External validation####

## Prepare: Mime1 expects consistent dataframes and a gene list. You’ve already got log2 expression, factors set.
## Set options: Use univariate Cox filter inside Mime1; include clinical covariates if supported as features.
# 1) Define your candidate genes (27 provided)
genelist <- c("ALB","APC","AQP1","ARF1","ATP7A","ATP7B","BACE1","CCS","CDK1","CDKN2A",
              "CYP1A1","DAXX","DLD","DLST","LIPT1","MTF1","NFAT5","NFE2L2","PAM","PARK7",
              "PDE3B","PPP1R15A","PRNP","SCO1","SCO2","SNCA","SNCG")

# 2) Rename and coerce survival columns; create ID if missing


library(caret)  # for stratified split

df_ex_vali <- read.table("merged_data_TCGA.txt", header=TRUE, sep="\t", check.names=FALSE)


normalize_two <- function(d){
  d <- as.data.frame(d)
  
  # Create ID if not present
  if(!"ID" %in% colnames(d)){
    if (!is.null(rownames(d)) && all(rownames(d) != "")) {
      d$ID <- rownames(d)
    } else if ("Patient" %in% colnames(d)) {
      d$ID <- d$Patient
    } else {
      d$ID <- seq_len(nrow(d))
    }
  }
  
  # Map futime/fustat to required names if needed
  if("OS.time" %in% colnames(d) == FALSE && "futime" %in% colnames(d)) {
    names(d)[names(d) == "futime"] <- "OS.time"
  }
  if("OS" %in% colnames(d) == FALSE && "fustat" %in% colnames(d)) {
    names(d)[names(d) == "fustat"] <- "OS"
  }
  
  # Type coercions
  d$OS.time <- as.numeric(d$OS.time)
  d$OS      <- as.numeric(d$OS)
  d$OS[is.na(d$OS)] <- 0
  d$OS[d$OS > 1] <- 1  # enforce binary
  
  # Reorder columns: ID, OS.time, OS, then the rest
  first3 <- c("ID","OS.time","OS")
  rest   <- setdiff(colnames(d), first3)
  d <- d[, c(first3, rest)]
  
  # Remove rows with invalid time
  d <- d[is.finite(d$OS.time) & d$OS.time > 0, ]
  
  return(d)
}

####----------------------1.C-index---------------------####
## 1.1绘制训练集和验证集在所有模型的C-index热图和条形图
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order =names(list_train_vali_Data),
               width = 0.35)
ggsave('./breast/1.1_Heatmap_C_index.png', width = 6.3, height = 12)

####---------------------------2.生存KM曲线----------------------------####
##选择最佳模型（RSF）绘制KM生存曲线，显示模型能够显著区分高风险和低风险患者。
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "RSF + Ridge",dataset = names(list_train_vali_Data)[i],
                              median.line = "hv", ## 中位生存线类型
                              cutoff = 0.5, ## 风险分组的截断值（中位）
                              conf.int = T, ## 显示置信区间
                              xlab="Day", ## X轴标签名
                              pval.coord=c(1000,0.9) ## p值的位置
  )
  )
}
aplot::plot_list(gglist=survplot,ncol=2)
ggsave('./Lung-png/2.1_survival_curve.png', width = 12.36, height = 5.77)


####------------------3.生存ROC曲线-------------------####
##绘制1年、3年和5年的ROC曲线，并比较最佳模型在不同时间点的预测性能
## 1年
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data, ## 训练集和验证集
                             mode = 'all',## 需与建模时一致，选择所用 
                             AUC_time = 1, ## 1年生存时间
                             auc_cal_method="KM" ## 使用Kaplan-Meier方法
)
## 3.1绘制1年AUC结果
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1], ## 选择验证集
            order= names(list_train_vali_Data), ## 排序
            width = 0.35, ## 热图单元格宽度
            year=1 ## 生存时间点
)
ggsave('./Lung-png/3.1_1year_AUC.png', width = 7.2, height = 12)
# Archive the entire result object for reproducibility
saveRDS(all.auc.1y, "all.auc.1y.rds")
#all.auc.3y<-readRDS("all.auc.3y.rds")
