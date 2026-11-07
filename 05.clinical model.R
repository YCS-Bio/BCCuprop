##########临床预测模型##########
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

setwd("D:\\Prog")  
rt=read.table("uniSigExpTime.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$futime[rt$futime<=0]=0.003

#对数据进行分组，构建模型
n=1     #分组的数目
for(i in 1:n){
	#############对数据进行分组############
	inTrain<-createDataPartition(y=rt[,2], p=0.5, list=F)
	train<-rt[inTrain,]
	test<-rt[-inTrain,]
	trainOut=cbind(id=row.names(train),train)
	testOut=cbind(id=row.names(test),test)
	
	#lasso回归
	x=as.matrix(train[,c(3:ncol(train))])
	y=data.matrix(Surv(train$futime,train$fustat))
	fit <- glmnet(x, y, family = "cox", maxit = 1000)
	cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
	coef <- coef(fit, s = cvfit$lambda.min)
	index <- which(coef != 0)
	actCoef <- coef[index]
	lassoGene=row.names(coef)[index]
	lassoSigExp=train[,c("futime", "fustat", lassoGene)]
	lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
	geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
	if(nrow(geneCoef)<2){next}

	#############构建COX模型#############
	multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
	multiCox=step(multiCox,direction = "both")
	multiCoxSum=summary(multiCox)

	#输出模型的公式
	outMultiTab=data.frame()
	outMultiTab=cbind(
		               coef=multiCoxSum$coefficients[,"coef"],
		               HR=multiCoxSum$conf.int[,"exp(coef)"],
		               HR.95L=multiCoxSum$conf.int[,"lower .95"],
		               HR.95H=multiCoxSum$conf.int[,"upper .95"],
		               pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
	outMultiTab=outMultiTab[,1:2]
	
	#输出train组风险文件
	riskScore=predict(multiCox,type="risk",newdata=train)
	coxGene=rownames(multiCoxSum$coefficients)
	coxGene=gsub("`","",coxGene)
	outCol=c("futime","fustat",coxGene)
	medianTrainRisk=median(riskScore)
	risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
	trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
		
	#输出test组风险文件
	riskScoreTest=predict(multiCox,type="risk",newdata=test)
	riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
	testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))
	
	#生存差异pvalue	
	diff=survdiff(Surv(futime, fustat) ~risk,data = train)
	pValue=1-pchisq(diff$chisq, df=1)
	diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
	pValueTest=1-pchisq(diffTest$chisq, df=1)

	#ROC曲线下的面积
	predictTime=1    #预测时间
	roc=timeROC(T=train$futime, delta=train$fustat,
	            marker=riskScore, cause=1,
	            times=c(predictTime), ROC=TRUE)
	rocTest=timeROC(T=test$futime, delta=test$fustat,
	            marker=riskScoreTest, cause=1,
	            times=c(predictTime), ROC=TRUE)	
	
	if((pValue<0.01) & (roc$AUC[2]>0.65) & (pValueTest<0.05) & (rocTest$AUC[2]>0.63)){
		#输出分组结果
		write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
		write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)
	  #lasso结果
	  write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
		pdf("lasso.lambda.pdf")
		plot(fit, xvar = "lambda", label = TRUE)
		dev.off()
		pdf("lasso.cvfit.pdf")
		plot(cvfit)
		abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
		dev.off()
	  #输出多因素结果
		write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
		write.table(trainRiskOut,file="risk.train.txt",sep="\t",quote=F,row.names=F)
		write.table(testRiskOut,file="risk.test.txt",sep="\t",quote=F,row.names=F)
		#所有样品的风险值
		allRiskOut=rbind(trainRiskOut, testRiskOut)
		write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)
		break
	}
}


##########桑基图##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(ggalluvial)
library(ggplot2)
library(dplyr)

CuproptosisCluFile="CuproptosisCluster.txt"
geneCluFile="geneCluster.txt"
riskFile="risk.all.txt" 

#读取输入文件
CuproptosisClu=read.table(CuproptosisCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
twoCluster=cbind(CuproptosisClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(risk))
rt=cbind(risk[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

#准备桑基图输入文件
rt=rt[,c("Cuproptosiscluster", "geneCluster", "risk", "fustat")]
colnames(rt)=c("Cuproptosiscluster", "geneCluster", "Risk", "Fustat")
rt[,"Fustat"]=ifelse(rt[,"Fustat"]==0, "Alive", "Dead")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#绘制桑基图
pdf(file="ggalluvial.pdf", width=6, height=5.5)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  #绘制线条颜色，forward说明线条颜色与前面柱状图一致，backward说明线条颜色与后面的柱状图一致
  geom_flow(width = 2/10,aes.flow = "forward") + 
  geom_stratum(alpha = .9,width = 2/10) +
  scale_fill_manual(values = mycol) +
  #size=3代表字体大小
  geom_text(stat = "stratum", size = 3,color="black") +
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  ggtitle("") + guides(fill = FALSE)                            
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

#绘制生存曲线函数
bioSurvival=function(inputFile=null, outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组生存差异，得到生存差异的p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("red", "blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 10,height =5.5)
  print(surPlot)
  dev.off()
}

#调用函数，绘制生存曲线
bioSurvival(inputFile="risk.train.txt", outFile="surv.train.pdf")
bioSurvival(inputFile="risk.test.txt", outFile="surv.test.pdf")
bioSurvival(inputFile="risk.all.txt", outFile="surv.all.pdf")


##########不同分型之间风险打分##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(limma)
library(ggpubr)
CuproptosisCluFile="Cuproptosiscluster.txt"      
geneCluFile="geneCluster.txt"    
scoreFile="risk.all.txt"       

#读取输入文件
CuproptosisClu=read.table(CuproptosisCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score$riskScore[score$riskScore>quantile(score$riskScore,0.99)]=quantile(score$riskScore,0.99)

#合并数据
twoCluster=cbind(CuproptosisClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

##失巢凋亡分型的箱线图
#设置比较组
data$Cuproptosiscluster=factor(data$Cuproptosiscluster, levels=levels(factor(data$Cuproptosiscluster)))
group=levels(factor(data$Cuproptosiscluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义图形颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$Cuproptosiscluster)))]

#绘制boxplot
boxplot=ggboxplot(data, x="Cuproptosiscluster", y="riskScore", color="Cuproptosiscluster",
                  xlab="Cuproptosiscluster",
                  ylab="Risk score",
                  legend.title="Cuproptosiscluster",
                  palette=bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图片
pdf(file="Cuproptosiscluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()

##基因分型的箱线图
#设置比较组
data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义图形颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]

#绘制boxplot
boxplot=ggboxplot(data, x="geneCluster", y="riskScore", color="geneCluster",
                  xlab="geneCluster",
                  ylab="Risk score",
                  legend.title="geneCluster",
                  palette=bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图片
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()


##########独立预后分析##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(survival)

#绘制森林图函数
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrLow[hrLow<0.001]=0.001
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width=6.5, height=4.8)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  LOGindex=2 
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,LOGindex^a1)
  dev.off()
}

#定义独立预后分析函数
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #数据合并
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
  #单因素独立预后分析
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")
  
  #多因素独立预后分析
  uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
  rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}

#调用函数，进行独立预后分析"
indep(riskFile = "risk.all.txt",
      cliFile = "clinical1.txt",
      uniOutFile = "uniCox.txt",
      multiOutFile = "multiCox.txt",
      uniForest = "uniForest.pdf",
      multiForest = "multiForest.pdf")


##########Nomogram、ROC、校准曲线##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(survival)
library(survminer)
library(timeROC)
library(rms)
library(regplot)
riskFile="risk.all.txt"    
cliFile="clinical2.txt"       

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#绘制列线图
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"), #图片展示的形式
             clickable=F, 
             title="",                      #图形的标题
             points=TRUE,                   #是否展示点
             droplines=TRUE,                #是否展示线段
             observation=rt[1,],            #展示病人的编号
             rank="sd",
             failtime = c(1,3,5),           #预测年限
             prfail = T)


#绘制ROC曲线
nomoRisk=predict(res.cox, data=rt, type="risk")
rt$nomoRisk=nomoRisk
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)
pdf(file="ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col='red',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='yellow',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='green',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("red","yellow","green"),lwd=2,bty = 'n')
dev.off()

#校准曲线
pdf(file="calibration.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ nomoRisk, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ nomoRisk, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5年校准曲线 
f <- cph(Surv(futime, fustat) ~ nomoRisk, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()


##########临床相关性ROC曲线##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(survivalROC)     
riskFile="risk.all.txt"          
cliFile="clinical1.txt"       

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime","fustat","riskScore")]

#读取临床数据文件
cli=read.table(cliFile,sep="\t",header=T,check.names=F,row.names=1)

#合并文件
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#绘制risk打分的ROC曲线
rocCol=rainbow(ncol(rt)-2)
aucText=c()
pdf(file="cliROC.pdf", width=6, height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=risk$futime, status=risk$fustat, marker=risk$riskScore, predict.time=1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制其他临床性状的ROC曲线
j=1
for(i in colnames(rt[,4:ncol(rt)])){
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
  j=j+1
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
  aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
}
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol)
dev.off()


##########免疫检查点分子##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(limma)
library(corrplot)

expFile="merge.txt"        
geneFile="gene.txt"         
riskFile="risk.all.txt"     

#读取基因表达文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
colnames(rt)=gsub("(.*?)\\_(.*?)", "\\2", colnames(rt))
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件，提取免疫检查点相关基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level=0.95)  #显著性标记的矩阵 

#绘制相关性图形
pdf(file="geneCor.pdf", width=7, height=7)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=1, pch=T,
         tl.col="black",
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


##########免疫功能##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(limma)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(reshape2)

expFile="merge.txt"       
gmtFile="immune.gmt"        
riskFile="risk.all.txt"    

#读取基因表达文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
colnames(rt)=gsub("(.*?)\\_(.*?)", "\\2", colnames(rt))
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

#读取数据集文件
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA 打分矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA 打分进行矫正
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)

#读取风险文件，获取高低风险
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=row.names(risk[risk$risk=="low",])
highSample=row.names(risk[risk$risk=="high",])
lowData=data[,lowSample]
highData=data[,highSample]
data=cbind(lowData, highData)
conNum=ncol(lowData)       #低风险样品
treatNum=ncol(highData)    #高风险样品
sampleType=c(rep(1,conNum), rep(2,treatNum))

#基因差异分析
sigVec=c()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ sampleType)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(i, Sig))
}
row.names(data)=sigVec

#热图可视化
Type=c(rep("Low risk",conNum), rep("High risk",treatNum))
Type=factor(Type, levels=c("Low risk", "High risk"))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf("heatmap.pdf", width=8, height=4.6)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()


##########免疫逃逸IPS##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(ggpubr)                    
tciaFile="TCIA.txt"                
scoreFile="risk.all.txt"   

#读取免疫治疗打分文件
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#读取失巢凋亡打分分组文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(ips), row.names(score))
ips=ips[sameSample, , drop=F]
score=score[sameSample, "risk", drop=F]
data=cbind(ips, score)

#设置比较组
data$group=factor(data$risk, levels=c("low", "high"))
group=levels(factor(data$risk))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对免疫治疗打分进行循环，分别绘制小提琴图
for(i in colnames(data)[1:(ncol(data)-1)]){
  rt=data[,c(i, "risk")]
  colnames(rt)=c("IPS", "risk")
  gg1=ggviolin(rt, x="risk", y="IPS", fill = "risk", 
               xlab="riskscore", ylab=i,
               legend.title="riskscore",
               palette=c("#0066FF", "#FF0000"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0(i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}


##########微卫星不稳定MSI##########
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())

library(plyr)
library(ggplot2)
library(ggpubr)
riskFile="risk.all.txt"     
cliFile="MSI.txt"          

#读取输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(risk), row.names(cli))
rt=cbind(risk[sameSample,,drop=F], cli[sameSample,,drop=F])
rt$MSI=factor(rt$MSI, levels=c("MSS", "MSI-L", "MSI-H"))
rt$risk=factor(rt$risk, levels=c("low", "high"))

#定义图形的颜色
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,"MSI"]))]

#统计高低风险组病人数目
rt1=rt[,c("MSI", "risk")]
df=as.data.frame(table(rt1))
#计算高低风险组的百分率
df=ddply(df, .(risk), transform, percent = Freq/sum(Freq) * 100)
#百分比位置
df=ddply(df, .(risk), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
#绘制百分比柱状图
p=ggplot(df, aes(x = factor(risk), y = percent, fill = MSI)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("riskScore")+ ylab("Percent weight")+  guides(fill=guide_legend(title="MSI"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()

#设置比较组
rt2=rt[,c("MSI", "riskScore")]
type=levels(factor(rt2[,"MSI"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#绘制箱线图
boxplot=ggboxplot(rt2, x="MSI", y="riskScore", fill="MSI",
                  xlab="",
                  ylab="Risk score",
                  legend.title="MSI",
                  palette=bioCol
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()
