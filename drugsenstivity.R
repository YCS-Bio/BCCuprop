install.packages(c("glmnet", "data.table"))
install.packages("BiocManager")
BiocManager::install("sva")
install.packages("devtools")
devtools::install_github("HuangLabUMN/oncoPredict")
sessionInfo()

library(oncoPredict)
library(dplyr)
library(tidyverse)
library(vroom)
library(reshape2)
library(ggpubr)

rm(list=ls())
options(stringsAsFactors=F, warn=-1)
groupFile="D://IC50/CuproptosisCluster2.txt"

exp <- read.table(
  "D:/IC50/DEGexpr.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

rownames(exp) <- trimws(rownames(exp))
exp <- exp[!duplicated(rownames(exp)), ]

exp <- as.matrix(exp)
mode(exp) <- "numeric"


library(httr)

CTRP2_Expr= readRDS(file="D://IC50/CTRP2_Expr (TPM, not log transformed).rds")
CTRP2_RES=readRDS(file="D:/IC50/CTRP2_Res.rds")
CTRP2_RES=exp(CTRP2_RES)
# write CSV with rownames (gene names) as first column
write.csv(CTRP2_Expr, file = "D:/IC50/CTRP2_Expr.csv", row.names = TRUE, quote = FALSE)
library(oncoPredict)

# 假设 CTRP2_RES 是 drugs x celllines 或 celllines x drugs，先查看结构
str(CTRP2_RES)
dim(CTRP2_RES)
head(rownames(CTRP2_RES))
head(colnames(CTRP2_RES))

# 统一为 drugs x samples（如果是反向请转置）
# 统计每个药物的非 NA 数量（按行）
drug_nonNA_counts <- rowSums(!is.na(as.matrix(CTRP2_RES)))
summary(drug_nonNA_counts)
table(cut(drug_nonNA_counts, breaks = c(0,5,10,20,50,Inf)))
# 查看某些药物的样本数
head(sort(drug_nonNA_counts, decreasing = FALSE))




pred <- calcPhenotype(
  trainingExprData = CTRP2_Expr,
  trainingPtype = CTRP2_RES,
  testExprData = exp,
  batchCorrect = "eb",
  powerTransformPhenotype = TRUE, 
  minNumSamples = 12,
  printOutput = TRUE,
  removeLowVaringGenesFrom = "homogenizeData"
)

head(pred)


#读入药敏文件
senstivity=read.csv("D:/IC50/calcPhenotype_Output/DrugPredictions.csv", header = T, sep = ",", check.names = F, row.names = 1)
group=read.table("D:/IC50/CuproptosisCluster2.txt", header = T, sep = "\t", check.names = F, row.names = 1)
  
#合并
sameSample=intersect(row.names(group),row.names(senstivity))
group=group[sameSample,"Cuproptosiscluster", drop=F]
senstivity=senstivity[sameSample,,drop=F]
senstivity[is.na(senstivity)]=0
senstivity=log2(senstivity+1)
df=senstivity
df$col=as.numeric(as.character(df$col))
options(scipen = 999)   # large positive value biases against scientific notation
# fixed decimal with 2 places for display
print(format(df$col, nsmall = 2, trim = TRUE))

# or use formatC for more control
print(formatC(df$col, format = "f", digits = 2))

# write to CSV with formatted strings (keeps numeric as text in file)
df$col_formatted <- formatC(df$col, format = "f", digits = 2)
write.csv(df, "out.csv", row.names = FALSE)
df=read.csv("D:/IC50/out.csv",header = T, sep = ",", check.names = F)

rtrt=cbind(group,df)
#设置比较组
rtrt$Cuproptosiscluster=factor(rtrt$Cuproptosiscluster, levels=c("1","2"))
type=levels(factor(rtrt[,"Cuproptosiscluster"]))
comp=combn(type,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#提取显著差异的药物
sigGene=c()
for(i in colnames(rtrt)[2:(ncol(rtrt))]){
  if(sd(rtrt[,2])<0.05){next}
  wilcoxTest=wilcox.test(rtrt[,i]~rtrt[,"Cuproptosiscluster"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.001){
    sigGene=c(sigGene,i)
  }
}
sigGene=c(sigGene,"Cuproptosiscluster")
rtrt=rtrt[,sigGene]

#把数据转换成ggplot2输入文件
rtrt=melt(rtrt,id.vars=c("Cuproptosiscluster"))
colnames(rtrt)=c("Cuproptosiscluster", "Gene", "Expression")



#设置比较组
#设置比较组
rtrt$Cuproptosiscluster=factor(rtrt$Cuproptosiscluster, levels=c("1","2"))
type=levels(factor(rtrt[,"Cuproptosiscluster"]))
comp=combn(type,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
boxplot=ggboxplot(rtrt, x="Gene", y="Expression", fill="Cuproptosiscluster",
                  xlab = "",
                  ylab = "Drug Senstivity",
                  legend.title="Cuproptosiscluster",
                  width=0.8,
                  palette = c("DodgerBlue1","Firebrick2"))+
  rotate_x_text(50)+
  stat_compare_means(aes(group = Cuproptosiscluster),
                     method="wilcox.test",
                     symum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),
                                     symbols=c("***","**","*","ns")),label = "p.signif")+
  theme(axis.text = element_text(face="bold.italic", colour = "#441718",size=16),
        axis.title = element_text(face="bold.italic", colour = "#441718",size=16),
        axis.line = element_blank(),
        plot.title = element_text(face="bold.italic", colour = "#441718",size=16),
        legend.text =element_text(face="bold.italic"),
        panel.border = element_rect(fill = NA, color="#35A79D",size = 1.5, linetype = "solid"),
        panel.background =  element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color="#CFD3D6",size=.5, linetype = "dotdash"),
        legend.title = element_text(face = "bold.italic", size = 13)
  )
#输出图片
pdf(file="drugSensitivity.pdf",width = 20, height = 8)
print(boxplot)
dev.off()




