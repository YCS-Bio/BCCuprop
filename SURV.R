#=============================================================
# Survival Analysis
#============================================================
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
setwd("D:/BioInf/expl")
workDir <- "D:/BioInf/expl"
library(survival)
library(survminer)
clusterFile="prediction.txt"
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
cat("Unique Cluster values (excluding NA):\n")
print(unique(rt$Cluster[!is.na(rt$Cluster)]))
cat("Table of Cluster values (including NA if any):\n")
print(table(rt$Cluster, useNA = "ifany"))
cat("Number of NAs in Cluster: ", sum(is.na(rt$Cluster)), "\n")
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
pdf(file="survival_BC.pdf",onefile = FALSE,width=12,height=5.5)
print(surPlot)
dev.off()