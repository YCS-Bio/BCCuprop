#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
setwd("D:/ML_REO/split0221")
results=readRDS("ConsensusClusterPlus_results2.rds")
consensus_mat <- results[[4]]$consensusMatrix

##=============================================================
# 2. 获取样本名（通常 consensusMatrix 的行名/列名就是样本名）
##=============================================================
samples <- names(results[[4]]$consensusClass)
if (is.null(samples)) samples <- colnames(consensus_mat)  # 有时行名为空，列名存在
colnames(consensus_mat)=samples
rownames(consensus_mat)=samples
# 3. 定义小簇的两个样本（请替换为实际样本名）
small_samples <- c("TCGA-E2-A1IN-01A-11R-A13Q-07", "TCGA-AN-A0AJ-01A-11R-A00Z-07")
# 4. 确保样本名存在于共识矩阵中
if (!all(small_samples %in% samples)) stop("样本名不在共识矩阵中")
# 5. 提取这两个样本的共分频率
# 共识矩阵中，位置 (i,j) 表示样本 i 和 j 被分到同一簇的频率。
# 我们要看这两个样本之间以及与所有其他样本的共分频率。

# 两个样本之间的共分频率
freq_between <- consensus_mat[small_samples[1], small_samples[2]]
cat("两个样本之间的共分频率：", freq_between, "\n")

# 两个样本与其他每个样本的平均共分频率（可选）
# 例如，样本1与其他所有样本的平均共分频率
others <- setdiff(samples, small_samples)
mean_freq1 <- mean(consensus_mat[small_samples[1], others])
cat("样本1与其他样本的平均共分频率：", mean_freq1, "\n")
mean_freq2 <- mean(consensus_mat[small_samples[2], others])
cat("样本2与其他样本的平均共分频率：", mean_freq2, "\n")
# 6. 可视化共识矩阵热图，高亮这两个样本
library(gplots)
# 创建一个颜色向量标记两个样本
col_colors <- rep("white", ncol(consensus_mat))
names(col_colors) <- samples
col_colors[small_samples] <- c("red", "blue")  # 分别标记不同颜色
# 绘制热图，添加列颜色条
heatmap.2(consensus_mat,
          trace = "none",
          main = paste("Consensus matrix k =", 4),
          ColSideColors = col_colors,
          RowSideColors = col_colors,
          key = TRUE,
          symkey = FALSE,
          density.info = "none")

dev.off()  # Close any open graphics device

##=========================
## 转置表达矩阵，使行为样本，列为基因
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
setwd("D:/ML_REO/split0221")
expr=read.table("expr_crg19_train_scaled.txt", header = T, row.names = 1, check.names = F)
expr_t <- t(expr)  # 现在行=样本，列=基因

# 计算 PCA
pca_res <- prcomp(expr_t, scale. = FALSE)  # 数据已标准化，无需再 scale
pca_df <- as.data.frame(pca_res$x[, 1:2])  # 取前两个主成分
pca_df$sample <- rownames(pca_df)
pca_df$cluster <- "other"  # 默认标记为其他
pca_df$cluster[rownames(pca_df) %in% small_samples] <- "small cluster"

# 绘图
library(ggplot2)
png("pca_plot.png", width = 800, height = 600)
ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("other" = "gray", "small cluster" = "red")) +
  theme_minimal() +
  ggtitle("PCA of samples, highlighting the two small-cluster samples")
dev.off()

##################
if (!require(Rtsne)) install.packages("Rtsne")
library(Rtsne)
set.seed(123)  # 固定随机种子
tsne_out <- Rtsne(expr_t, perplexity = 30, check_duplicates = FALSE)
tsne_df <- data.frame(tsne_out$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$sample <- rownames(expr_t)
tsne_df$cluster <- "other"
tsne_df$cluster[rownames(expr_t) %in% small_samples] <- "small cluster"

png("tsne_plot.png", width = 800, height = 600)
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("other" = "gray", "small cluster" = "red")) +
  theme_minimal() +
  ggtitle("t-SNE of samples, highlighting the two small-cluster samples")
dev.off()

##============================================================================
# 移除 并重新计算后续分析所需的所有对象
##============================================================================

# 假设这两个样本名如下（请根据实际替换）
small_samples <- c("TCGA-E2-A1IN-01A-11R-A13Q-07", 
                   "TCGA-AN-A0AJ-01A-11R-A00Z-07")

# 确认它们在训练集中
# expr_tcga 是原始表达矩阵（基因×样本），或者 expr_scaled 是标准化后的矩阵
# 这里以 expr_scaled 为例
if (!all(small_samples %in% colnames(expr_scaled))) {
  stop("样本名不在表达矩阵中，请检查")
}

# 原始标准化后的矩阵（基因×样本）
# 保留不在 small_samples 中的样本
expr_scaled=read.table("expr_crg19_train_scaled.txt", header = T, row.names = 1, check.names = F)
expr_scaled_clean <- expr_scaled[, !colnames(expr_scaled) %in% small_samples, drop = FALSE]

# 检查新矩阵的维度
dim(expr_scaled_clean)  # 应为 基因数 × (838-2) = 基因数 × 836

# 移除这两个样本的标签
clusters_clean=read.table("CuproptosisCluster4_1.txt", header = T, row.names=1, check.names = F)

# 检查剩余的簇标签分布
table(clusters_clean)
# 应该只剩下 1,2,3 三个簇（簇4已消失）

##重新计算标准化参数（scaler）——基于清理后的训练集
# 对清理后的表达矩阵重新标准化
# 注意：如果 expr_scaled 已经是标准化后的，需要先“反标准化”回原始值？不，直接基于原始表达矩阵重新标准化更安全。
# 假设您有原始表达矩阵 expr_raw（基因×样本，未标准化）

# 方法A：如果您有原始表达矩阵 expr_raw
expr_raw_no_log=read.table("expr_crg19_train.txt",header = T, row.names = 1,check.names = F)
expr_raw=log2(expr_raw_no_log+1)
expr_raw_clean <- expr_raw[, !colnames(expr_raw) %in% small_samples, drop = FALSE]
scaler_new <- scale(t(expr_raw_clean))
center_new <- attr(scaler_new, "scaled:center")
scale_new <- attr(scaler_new, "scaled:scale")
expr_scaled_new <- scaler_new  # 新的标准化矩阵

# 方法B：如果您没有原始表达矩阵，只有 expr_scaled（已标准化）
# 这种情况下无法准确反标准化，但可以假设移除两个样本对整体分布影响极小，直接使用原来的 scaler？
# 不推荐！必须基于新数据重新计算。
# 因此，请务必使用原始表达矩阵重新标准化。
##=============================
# 5. 重新计算三个大簇的质心
##=============================
# 获取清理后的样本对应的新聚类标签（确保顺序一致）
# 假设 clusters_clean 的顺序与 expr_scaled_clean 的列名一致
common_samples <- intersect(colnames(expr_scaled_clean), rownames(clusters_clean))
expr_scaled_clean <- expr_scaled_clean[, common_samples, drop = FALSE]
names(clusters_clean) <- colnames(expr_scaled_clean)

# 计算每个簇的质心（行=基因，列=簇）
new_centroids <- sapply(sort(unique(clusters_clean)), function(k) {
  samples_k <- names(clusters_clean[clusters_clean == k])
  rowMeans(expr_scaled_clean[, samples_k, drop = FALSE])
})
# 设置列名为簇标签（例如 "1","2","3"）
colnames(new_centroids) <- sort(unique(clusters_clean))
##=====================================
# 6. 保存清理后的所有对象
##=====================================
# 保存新的标准化参数
saveRDS(list(center = center_new, scale = scale_new), file = "scaler_clean.rds")

# 保存新的质心
saveRDS(new_centroids, file = "centroids_clean.rds")

# 保存新的聚类标签（可选）
write.table(data.frame(sample = names(clusters_clean), 
                       cluster = clusters_clean),
            file = "clusters_clean.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("清理完成！新训练集样本数：", ncol(expr_scaled_clean), "\n")
cat("簇大小：\n")
print(table(clusters_clean))


##==================================
# test
##==================================
# 加载清理后的对象
scaler_clean <- readRDS("scaler_clean.rds")
centroids_clean <- readRDS("centroids_clean.rds")

# 标准化测试集（假设 test_expr 是原始表达矩阵，基因×样本）
test=read.table("expr_crg19_test.txt", header = T, row.names = 1,check.names = F)
test_expr=log2(test+1)
test_expr_t=t(test_expr)
test_scaled <- scale(test_expr_t, 
                     center = scaler_clean$center[colnames(test_expr_t)], 
                     scale = scaler_clean$scale[colnames(test_expr_t)])

# 分配测试集样本到三个大簇
# 假设 centroids_clean 的列是簇名，行是特征
test_clusters <- apply(test_scaled, 1, function(sample_vec) {
  # 计算该样本与每个簇中心的相关系数
  sims <- apply(centroids_clean, 2, function(centroid_vec) {
    cor(sample_vec, centroid_vec, method = "pearson")
  })
  # 返回相关系数最大的簇名
  colnames(centroids_clean)[which.max(sims)]
})
# 保存 test_clusters 到 CSV
# 如果 test_clusters 是向量
write.csv(
  data.frame(Sample = rownames(test_scaled), Cluster = test_clusters),
  file = "test_clusters.csv",
  row.names = FALSE
)
table(test_clusters)
#=============================================================
# Survival Analysis
#============================================================
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
setwd("D:/ML_REO/split0221")
workDir <- 'D:/ML_REO/split0221'
library(survival)
library(survminer)
clusterFile="test_clusters.txt"
cliFile="time.txt"
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
#colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
#生存差异统计
# --- DEBUGGING START ---
cat("Number of rows in rt: ", nrow(rt), "\n")
cat("Unique Cuproptosiscluster values (excluding NA):\n")
print(unique(rt$Cuproptosiscluster[!is.na(rt$Cuproptosiscluster)]))
cat("Table of Cuproptosiscluster values (including NA if any):\n")
print(table(rt$Cuproptosiscluster, useNA = "ifany"))
cat("Number of NAs in Cuproptosiscluster: ", sum(is.na(rt$Cuproptosiscluster)), "\n")
cat("Unique fustat values (excluding NA):\n")
print(unique(rt$fustat[!is.na(rt$fustat)]))
cat("Table of fustat values (including NA if any):\n")
print(table(rt$fustat, useNA = "ifany"))
cat("Number of NAs in fustat: ", sum(is.na(rt$fustat)), "\n")
# --- DEBUGGING END ---


length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))
#绘制生存曲线
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit,
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Cluster",
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   xlim = c(0,14),
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival_test.pdf",onefile = FALSE,width=12,height=5.5)
print(surPlot)
dev.off()

#=============================================================================
#GSVA
#=============================================================================
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
setwd("D:/ML_REO/split0221")
workDir <- 'D:/ML_REO/split0221' 
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
expFile="TCGA_TPM.txt"             
clusterFile="test_clusters.txt"     
gmtFile="c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt"         

#读取表达输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", row.names = 1, check.names = F)
rt=as.matrix(rt)
# rownames(rt)=rt[,1]
exp=rt
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=log2(data+1)

#GSVA富集分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvapar=gsvaParam(data, 
                  geneSets, 
                  kcdf='Gaussian',
                  absRanking=FALSE)
gsvaResult <- gsva(gsvapar)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)
gsvaOut=read.table("gsvaOut.txt", sep="\t", header = T,check.names=F)
gsvaResult=gsvaOut
#读取cluster文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
gsvaResult1=gsvaResult

#数据合并
gsvaResult=t(gsvaResult)
#rownames(gsvaResult)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
colnames(gsvaResult) <- gsvaResult[1, ]     # 第一行作为列名
gsvaResult <- gsvaResult[-1, ]               # 删除第一行
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
gsvaCluster=cbind(gsvaCluster, Project)

#差异分析
adj.P.Val.Filter=1e-8
allType=as.vector(gsvaCluster$Cluster)
comp=combn(levels(factor(allType)), 2)

for(i in 1:ncol(comp)){
  #样品分组
  treat=gsvaCluster[gsvaCluster$Cluster==comp[2,i],]
  con=gsvaCluster[gsvaCluster$Cluster==comp[1,i],]
  data=rbind(con, treat)
  
  #expr <- data[, !(names(data) %in% "Cluster")]
  #expr <- expr[, !(names(data) %in% "Project")]
  #expr <- as.matrix(expr)
  #if (!is.numeric(expr[1, 1])) {
  #  stop("Expression matrix contains non-numeric values.")
  #}
  
  #差异分析
  Type=as.vector(data$Cluster)
  ann = data.frame(Cluster = Type)
  rownames(ann) = rownames(data)
  data=data[,-c((ncol(data)-1), ncol(data))]
  data <- as.matrix(data)
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  
  fit=lmFit(data, design)
  
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #输出所有通路的差异情况
  allDiff=topTable(fit2, adjust='fdr', number=200000)
  allDiffOut=rbind(id=colnames(allDiff), allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #输出显著的差异 (按|logFC|排序, adj.P.Val < 1e-8)
  diffSig=allDiff[with(allDiff, adj.P.Val < adj.P.Val.Filter), ]
  diffSig=diffSig[order(abs(diffSig$logFC), decreasing=TRUE), ]
  diffSigOut=rbind(id=colnames(diffSig), diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  
  #聚类颜色
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  CuproptosisCluCol=bioCol[1:length(levels(factor(allType)))]
  names(CuproptosisCluCol)=levels(factor(allType))
  ann_colors[["Cluster"]]=CuproptosisCluCol[c(comp[1,i], comp[2,i])]
  
  #绘制差异通路热图 (取|logFC|排序前20个)
  termNum=20
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  
  pdf(file=paste0(contrast,".heatmap1.pdf"), height=6, width=12)
  pheatmap(hmExp,
           annotation=ann,
           annotation_colors=ann_colors,
           color=colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
           cluster_cols=F,
           show_colnames=F,
           gaps_col=as.vector(cumsum(table(Type))),
           scale="row",
           fontsize=10,
           fontsize_row=7,
           fontsize_col=10)
  dev.off()
}


write.table(diffSigOut, file="diffSigout.txt", sep="\t", quote=F, col.names=F)

# 输出显著的差异 (按|logFC|排序, adj.P.Val < 1e-8)
diffSig = allDiff[with(allDiff, adj.P.Val < adj.P.Val.Filter), ]
diffSig = diffSig[order(abs(diffSig$logFC), decreasing=TRUE), ]

# 保存所有显著结果
diffSigOut = rbind(id=colnames(diffSig), diffSig)
write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)

# 取前20个通路并保存到CSV
top20 = head(diffSig, 20)
write.csv(top20, file=paste0(contrast, ".top20.csv"), row.names=TRUE)

# 绘制差异通路热图 (取|logFC|排序前20个)
termNum=20
diffTermName=as.vector(rownames(diffSig))
diffLength=length(diffTermName)
if(diffLength<termNum){termNum=diffLength}
hmGene=diffTermName[1:termNum]
hmExp=data[hmGene,]

# 将过长的通路名拆成两行 (在空格处插入换行符)
rownames(hmExp) <- sapply(rownames(hmExp), function(x) {
  # 按空格拆分，前半部分和后半部分用换行符连接
  parts <- strsplit(x, " ")[[1]]
  if(length(parts) > 1){
    mid <- ceiling(length(parts)/2)
    paste(c(paste(parts[1:mid], collapse=" "),
            paste(parts[(mid+1):length(parts)], collapse=" ")), 
          collapse="\n")
  } else {
    x
  }
})

# 保存为PNG
png(file=paste0(contrast,".heatmap.png"), height=6, width=12, units="in", res=300)
pheatmap(hmExp,
         annotation=ann,
         annotation_colors=ann_colors,
         color=colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
         cluster_cols=F,
         show_colnames=F,
         gaps_col=as.vector(cumsum(table(Type))),
         scale="row",
         fontsize=10,
         fontsize_row=7,
         fontsize_col=10)
dev.off()

