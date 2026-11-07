#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
# Install if needed
if (!requireNamespace("ConsensusClusterPlus", quietly=TRUE)) {
  BiocManager::install("ConsensusClusterPlus")
}
library(ConsensusClusterPlus)
setwd("D:/ML_REO/split0221")
# 1 Load data (columns = samples)
expr <- read.delim("expr_crg19_train.txt", check.names=FALSE)
gene_col <- colnames(expr)[1]
mat <- as.matrix(expr[ , -1])        # remove gene id column
rownames(mat) <- expr[[1]]

# 2 Preprocess
mat <- log2(mat + 1)                 # transform
# optional: filter top variable genes
vars <- apply(mat, 1, var)
mat <- mat[order(vars, decreasing=TRUE)[1:19], ]

#  z score by gene
mat_z <- t(scale(t(mat), center = TRUE, scale = TRUE))
write.table(mat_z, file="expr_crg19_train_scaled.txt", sep="\t", quote=F, col.names=T)

# save scaled:center & scale

center <- attr(mat_z, "scaled:center")
scale  <- attr(mat_z, "scaled:scale")
scaler <- list(center = center, scale = scale)
saveRDS(scaler, "scaler.rds")






# 3 Run ConsensusClusterPlus
results <- ConsensusClusterPlus(
  d = mat_z,
  maxK = 10,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  innerLinkage = "average",
  finalLinkage = "average",
  distance = "pearson",
  seed = 1234,
  title = "consensus_patient_1097",
  plot = "png",
  writeTable = TRUE
)
# Save the entire results object from ConsensusClusterPlus
saveRDS(results, file="ConsensusClusterPlus_results2.rds")
cat("ConsensusClusterPlus_results2.rds saved successfully.\n")

#输出分型结果
clusterNum=9     #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cuproptosiscluster")
letter=c("A","B","C","D","E","F","G","H","I","J")
uniqClu=levels(factor(cluster$Cuproptosiscluster))
cluster$Cuproptosiscluster=letter[match(cluster$Cuproptosiscluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="CuproptosisCluster9.txt", sep="\t", quote=F, col.names=F)

##=============================================
#  热图
##=============================================

#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
library(pheatmap)
setwd("D:/ML_REO/split0221")
workDir <- 'D:/ML_REO/split0221'
expFile="expr_crg19_train.txt"
clusterFile="CuproptosisCluster9.txt"
cliFile="cli_986.txt"

#读取输入文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
# The 'exp_scaled_data.txt' already has samples as rows, so no need to transpose here.
exp=t(exp) # <--- REMOVED THIS LINE
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#合并表达和分型数据
sameSample=intersect(row.names(exp), row.names(cluster))

expCluster=exp[sameSample,,drop=F]

#cli=read.csv(cliFile, header=T, check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample1=intersect(row.names(expCluster), row.names(cli))
expCluster1=expCluster[sameSample1,,drop=F]
cli1=cli[sameSample1,,drop=F]
data=cbind(expCluster1, cli1)

# 确保 data 是数据框而不是矩阵
data <- as.data.frame(data)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

# 'expCluster' from previous steps contains gene expression and Cuproptosiscluster
# 'cli' from previous steps contains clinical data

# Extract gene expression matrix separately
gene_expression_matrix <- expCluster

# Create the annotation data frame, combining Cuproptosiscluster from expCluster and cli
# Ensure samples are correctly aligned
common_samples = intersect(rownames(gene_expression_matrix), rownames(cli))

# Ensure Cuproptosiscluster is taken from expCluster and cli is properly aligned
annotation_df_cluster = data.frame(Cuproptosiscluster = cluster[common_samples, "Cuproptosiscluster", drop=FALSE])
rownames(annotation_df_cluster) = common_samples

annotation_df = cbind(annotation_df_cluster, cli[common_samples, , drop=FALSE])

# Order samples by Cuproptosiscluster
order_indices = order(annotation_df$Cuproptosiscluster)
annotation_df = annotation_df[order_indices, , drop=FALSE]

# Order the gene expression matrix columns (samples) to match the annotation order
mat_for_heatmap_untransposed = gene_expression_matrix[rownames(annotation_df), ] # Samples as rows, genes as columns
mat_for_heatmap = t(mat_for_heatmap_untransposed) # Transpose to genes as rows, samples as columns for pheatmap

# 'Type' will be the annotation data frame for pheatmap's annotation_col
Type = annotation_df

#聚类颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CuproptosisCluCol=bioCol[1:length(levels(factor(Type$CuproptosisCluster)))]
names(CuproptosisCluCol)=levels(factor(Type$CuproptosisCluster))
ann_colors[["CuproptosisCluster"]]=CuproptosisCluCol


#热图可视化
pheatmap(mat_for_heatmap, # Use the prepared numeric expression matrix
         annotation_col=Type, # Use the combined annotation data frame for columns
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6,
         filename="heatmap9.pdf", # Save directly using filename argument
         width=8, height=6) # Specify dimensions here

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
clusterFile="CuproptosisCluster7.txt"
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


length=length(levels(factor(rt$Cuproptosiscluster)))
diff=survdiff(Surv(futime, fustat) ~ Cuproptosiscluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cuproptosiscluster, data = rt)
#print(surv_median(fit))
#绘制生存曲线
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit,
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Cuproptosiscluster",
                   legend.labs=levels(factor(rt[,"Cuproptosiscluster"])),
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
pdf(file="survival7.pdf",onefile = FALSE,width=12,height=5.5)
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
clusterFile="CuproptosisCluster4.txt"     
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

#读取cluster文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
gsvaResult1=gsvaResult

#数据合并
gsvaResult=t(gsvaResult)
#rownames(gsvaResult)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
gsvaCluster=cbind(gsvaCluster, Project)

#差异分析
adj.P.Val.Filter=1e-8
allType=as.vector(gsvaCluster$Cuproptosiscluster)
comp=combn(levels(factor(allType)), 2)

for(i in 1:ncol(comp)){
  #样品分组
  treat=gsvaCluster[gsvaCluster$Cuproptosiscluster==comp[2,i],]
  con=gsvaCluster[gsvaCluster$Cuproptosiscluster==comp[1,i],]
  data=rbind(con, treat)
  
  #差异分析
  Type=as.vector(data$Cuproptosiscluster)
  ann = data.frame(Cuproptosiscluster = Type)
  rownames(ann) = rownames(data)
  data=t(data[,-c((ncol(data)-1), ncol(data))])
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
  ann_colors[["Cuproptosiscluster"]]=CuproptosisCluCol[c(comp[1,i], comp[2,i])]
  
  #绘制差异通路热图 (取|logFC|排序前20个)
  termNum=20
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  
  pdf(file=paste0(contrast,".heatmap.pdf"), height=6, width=12)
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

