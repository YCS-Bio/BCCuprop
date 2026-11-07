#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")


#BiocManager::install("limma")
library(S4Vectors)
packageVersion("S4Vectors")  # Should be >= 0.47.6 (e.g., 0.48.0)
library(Seqinfo)
packageVersion("Seqinfo")    # Should be 1.0.0 or higher

library(limma)           
expFile="merge.txt"      
geneFile="gene.txt"     
setwd("D:/tcga/07.AnoikisExp") 

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
write.table(out,file="AnoikisGeneExp.txt",sep="\t",quote=F,col.names=F)
