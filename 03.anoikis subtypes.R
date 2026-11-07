##########预后相关失巢凋亡基因表达矩阵##########
library(limma)           
expFile="merge.txt"      
geneFile="PrognosticAnoikis.txt"     
setwd("D:\\AnoikisCRC\\03.anoikis subtypes") 

#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#获取失巢凋亡基因的表达量
gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#输出结果
out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="PrognosticAnoikisExp.txt",sep="\t",quote=F,col.names=F)


##########失巢凋亡分型##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

workDir="D:\\AnoikisCRC\\03.anoikis subtypes"    
library(ConsensusClusterPlus)     
expFile="PrognosticAnoikisExp.txt"

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#聚类
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


#输出分型结果
clusterNum=3      #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cuproptosiscluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$Cuproptosiscluster))
cluster$Cuproptosiscluster=letter[match(cluster$Cuproptosiscluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="Cuproptosiscluster.txt", sep="\t", quote=F, col.names=F)


##########热图##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(pheatmap)        
expFile="PrognosticAnoikisExp.txt"         
clusterFile="Cuproptosiscluster.txt"     
cliFile="clinical.txt"

#读取输入文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#合并表达和分型数据
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
rownames(expCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expCluster))
expCluster=cbind(expCluster, Project)

#合并临床数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

#提取热图数据
data=data[order(data$Cuproptosiscluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
data=t(data[,1:ncol(exp)])

#聚类颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
AnoikisCluCol=bioCol[1:length(levels(factor(Type$Cuproptosiscluster)))]
names(AnoikisCluCol)=levels(factor(Type$Cuproptosiscluster))
ann_colors[["Cuproptosiscluster"]]=AnoikisCluCol

#热图可视化
pdf("heatmap.pdf", height=6, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()


##########生存分析##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(survival)
library(survminer)
clusterFile="Cuproptosiscluster.txt"     
cliFile="time.txt"               

#读取输入文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#生存差异统计
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
pdf(file="survival.pdf",onefile = FALSE,width=12,height=5.5)
print(surPlot)
dev.off()


##########GSVA富集分析###########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

setwd("D:\\1-7cuproptosisenrichment") 

library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
expFile="TCGA_TPM.txt"             
clusterFile="CuproptosisCluster2.txt"     
gmtFile="c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt"         

#读取表达输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]



dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=log2(data+0.1)

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

#数据合并
gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
gsvaCluster=cbind(gsvaCluster, Project)

#差异分析
adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$Cuproptosiscluster)
comp=combn(levels(factor(allType)), 2)
for(i in 1:ncol(comp)){
  #样品分组
  treat=gsvaCluster[gsvaCluster$Cuproptosiscluster==comp[2,i],]
  con=gsvaCluster[gsvaCluster$Cuproptosiscluster==comp[1,i],]
  data=rbind(con, treat)
  #差异分析
  Type=as.vector(data$Cuproptosiscluster)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #输出所有通路的差异情况
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #输出显著的差异
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  
  #聚类颜色
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  CuproptosisCluCol=bioCol[1:length(levels(factor(allType)))]
  names(CuproptosisCluCol)=levels(factor(allType))
  ann_colors[["Cuproptosiscluster"]]=CuproptosisCluCol[c(comp[1,i], comp[2,i])]
  
  #绘制差异通路热图
  termNum=20
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  pdf(file=paste0(contrast,".heatmap.pdf"),height=6,width=14)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
           cluster_cols =F,
           show_colnames = F,
           gaps_col=as.vector(cumsum(table(Type))),
           scale="row",
           fontsize = 10,
           fontsize_row=7,
           fontsize_col=7)
  dev.off()
}


##########GO富集分析##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library("ggpubr")
library(ComplexHeatmap)

pvalueFilter=0.05    
qvalueFilter=0.05   


colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

rt=read.table("PrognosticCup.txt", header=T, sep="\t", check.names=F)

#提取预后基因名称，将基因名臣转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#输出显著富集的结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#绘制分类柱状图
data=GO %>% group_by(ONTOLOGY) %>% slice_head(n=10)
pdf(file="barplot.color.pdf", width=8, height=6.5)
ggbarplot(data, x="Description", y="Count", fill = "ONTOLOGY", color = "white",
          xlab="Term",
          orientation = "horiz",   
          palette = "aaas",   
          legend = "right", 
          sort.val = "asc",   
          sort.by.groups=TRUE)+ 
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))
dev.off()


##########ssGSEA富集分析##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
expFile="TCGA_TPM.txt"           
gmtFile="immune.gmt"         
clusterFile="Cuproptosiscluster2.txt"       

#读取表达输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=log2(data+0.1)
#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
##old method:ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

gseapar=ssgseaParam(data,
                    geneSets
                    )
ssgseaScore=gsva(gseapar)


#对ssGSEA打分进行校正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("Cuproptosiscluster"))
colnames(data)=c("Cuproptosiscluster", "Immune", "Fraction")

#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"Cuproptosiscluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="Cuproptosiscluster", 
            ylab="Immune infiltration",
            xlab="",
            legend.title="Cuproptosiscluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot.pdf", width=8, height=6.5)    
p+stat_compare_means(aes(group=Cuproptosiscluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

