source("http://bioconductor.org/biocLite.R")
biocLite("limma")
install.packages("magrittr")
install.packages("ggplot2")
install.packages("corrplot")
install.packages("gplots")
install.packages("survminer")
install.packages("survival")
install.packages("glmnet")
#####Limma包进行癌和癌旁差异位点筛选
library(limma)
library(magrittr)
setwd("E:\\geo_learn\\statistics_analysis")
exp=read.csv(file = "exp_prepared.csv")
#exp[1:10,1:10]
rownames(exp)=exp[,"probe"]
#head(rownames(exp))
exp=exp[,-1]
#exp[1:10,1:10]
samps <- colnames(exp) 
#head(samps)
samps<-substr(samps,14,16)
#head(samps)
samps<-factor(samps,labels=c("N","T"))
#head(samps)
#build the design matrix for the lm
design <- model.matrix(~0+samps) 
#head(design)
colnames(design) <- c("Normal","Tumor") 
#head(design)
#bulid the contrast model
cont.matrix<-makeContrasts(Tumor-Normal,levels=design)

#bulid the lm model
fit <- lmFit(exp, design)
#compute the difference between the cancer and normal
fit2 <- contrasts.fit(fit, cont.matrix) 
#Empirical Bayes to calculate the statistic
fit2 <- eBayes(fit2)
# variable seletion
all_data <- topTable(fit2, coef=1, adjust.method="fdr", p.value=1, lfc=0, number=300000)
selected_data <-topTable(fit2, coef=1, adjust.method="fdr", p.value=0.0001, lfc=1, number=30000)
write.csv(selected_data,file = "selected_data.csv")
#####火山图
library(ggplot2)
#define a function of DrawVolcanoplot 
DrawVolcanoplot <- function(a,targetlogFC = 1, targetP.Value = 0.05, target.Points = NULL){
  targetPoints <- subset(a,a$X %in% target.Points)
  P.Value <- c(a$adj.P.Val)
  logFC <- c(a$logFC)
  df <- data.frame(P.Value, logFC)
  df.G <- subset(df, logFC < -targetlogFC& P.Value < targetP.Value) #define Green
  df.G <- cbind(df.G, rep("low", nrow(df.G)))
  colnames(df.G)[3] <- "Colors"
  df.B <- subset(df, (logFC >= -targetlogFC & logFC <= targetlogFC) | P.Value >= targetP.Value) #define Black
  df.B <- cbind(df.B, rep("not significant", nrow(df.B)))
  colnames(df.B)[3] <- "Colors"
  df.R <- subset(df, logFC > targetlogFC & P.Value < targetP.Value) #define Red
  df.R <- cbind(df.R, rep("high", nrow(df.R)))
  colnames(df.R)[3] <- "Colors"
  df.t <- rbind(df.G, df.B, df.R)
  df.t$Colors <- as.factor(df.t$Colors)
  p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value)) ) +
    geom_point(alpha = 0.2, size = 5,aes(color= Colors)) +
    
    xlim(c(-4, 4)) + ylim(c(0, 30)) +
    geom_vline(xintercept = c(-targetlogFC,targetlogFC),linetype=2) +
    geom_hline(yintercept = -log10(targetP.Value),linetype=2) +
    scale_color_manual(values = c("blue", "grey", rgb(185,0,0,maxColorValue=255))) + # R浣跨敤RGB棰滆壊锛岄粯璁ゆ槸0~1
    labs(x=expression(log[2](FC)), y=expression( -log[10](P.Value))) +
    theme(axis.title.x=element_text(size=20), 
          axis.text.x=element_text(size=15)) +
    theme(axis.title.y=element_text(size=20),
          axis.text.y=element_text(size=15)) +
    geom_text(data=targetPoints,aes(x = logFC, y = -log10(adj.P.Val)),label=targetPoints$X,check_overlap = F,hjust = 0, nudge_x = 0.5,color="black",size=8) +
    
    geom_point(data=targetPoints,aes(x = logFC, y = -log10(adj.P.Val)),color="black",size = 5)
  return(p)
}

#set the target log(FC);and target adj.P.value
DrawVolcanoplot(all_data,1,0.0001)

#####直方图
#select variable needed in selected_data
hist.data=selected_data[,c(1,4,5)]
#head(hist.data)
hist.data$y=-log10(hist.data$adj.P.Val)
#head(hist.data)
hist.data$pos=hist.data$logFC>=0
#head(hist.data)
hist.data$number=sample(c(1:84),84,replace = F)
#head(hist.data)
ggplot(hist.data, aes(x=number,y=y, fill=pos)) +geom_bar(stat="identity",position = position_dodge(width=0.1))+scale_fill_manual(values=c("blue", rgb(185,0,0,maxColorValue=255)))+
  theme(panel.grid.minor=element_blank(),
        axis.line=element_line(size=0.5))+labs(list(x = "n=84", y = "-Log10 of P"))+guides(color=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(size = 15,face = "bold"))+
  theme(legend.text= element_text(size=20, face= "bold"))+
  theme(legend.title= element_text(size=20, face= "bold"))+
  theme(axis.text.y = element_text(size = 15,face = "bold"))+
  theme(axis.title.x = element_text(size = 20, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(legend.position = "none")

#####相关系数图
#data frame transformed
#exp[1:10,1:10]
data.exp=as.data.frame(t(exp))
#data.exp[1:10,1:10]
#selected the data by limma
data_X_Y=data.exp[,c(rownames(selected_data))]
#head(rownames(selected_data))
#data_X_Y[1:10,1:10]
#special condition
i=1;while (i<dim(data_X_Y)[2]){
  data_X_Y[,i]=2^(as.numeric(as.matrix(data_X_Y[,i])))-1;
  i=i+1
}

library(corrplot)
#correlation coefficient
corr=cor(data_X_Y)
#absolute value
corr=abs(corr)
#plot
corrplot(corr = corr,method = "color",type="upper",order = "hclust",tl.pos="n",cl.lim=c(0,1))


#####LASSO logistic regression
library(glmnet)
#add the outcome
data_X_Y$object=c(rep(0,49),rep(1,49))
head(data_X_Y)
#matix x for the logistic lasso
x=data.matrix(data_X_Y[,-dim(data_X_Y)[2]])
#names(data_X_Y)
#dim(data_X_Y)[2]]
#matix y for the logistic lasso
y=as.matrix(data_X_Y[,dim(data_X_Y)[2]])
#logistic lasso function
fit1=glmnet(x,y,family="binomial")
plot(fit1)
#cross validation
set.seed(1234)
cv.fit1=cv.glmnet(x,y,family="binomial")
plot(cv.fit1); cv.fit1$lambda.min;cv.fit1$lambda.1se

#get coefficients
predict(cv.fit1,type='coefficients',s=cv.fit1$lambda.min)
#get the valiables selected by lasso logistic
a=predict(cv.fit1,type='nonzero',s=cv.fit1$lambda.min)
colnames(data_X_Y)[a$X1]
#plot the change of the coefficients
plot(cv.fit1$glmnet.fit,xvar = "lambda",label = T)
abline(v=log(c(cv.fit1$lambda.min,cv.fit1$lambda.1se)),lty=2)


#####热图
library(gplots)
#get the data for heatmap
data.heatmap=t(data_X_Y[,a$X1][,-length(names(data_X_Y))])
#head(data.heatmap)
#get the characteristics of the tumor or normal
strings=substring(colnames(data.heatmap),14,15)
#define the color of sample
color.map2<- function(strings) { if (strings=="11") "#FF0000" else "#0000FF" }
patientcolors2 <- unlist(lapply(strings, color.map2))
#plot
heatmap.2(data.heatmap, Colv = F,col=redgreen(75),scale="row", ColSideColors=patientcolors2,key=TRUE,symkey=FALSE, density.info="none", trace="none", dendrogram="none",cexRow=0.5)


#####LASSO cox regression
#input the data
exp_clinical=read.csv(file = "exp_clinical.csv")[,-1]
x=data.matrix(exp_clinical[,-c(1,2)])
y=cbind(time=exp_clinical$time,status=exp_clinical$event)
fit2=glmnet(x,y,family="cox")
plot(fit2)
set.seed(1234)
cv.fit2=cv.glmnet(x,y,family="cox")
plot(cv.fit2);cv.fit2$lambda.min; cv.fit2$lambda.1se
predict(cv.fit2,type='coefficients',s=cv.fit2$lambda.min)
plot(cv.fit2$glmnet.fit,xvar = "lambda")
abline(v=log(cv.fit2$lambda.min),lty=2)

######Kaplan Meier analysis
library(survminer);library(survival)
attach(exp_clinical)
exp_clinical$risk_score=X1*(1.662476e-03)+X2*(-1.621524e-03)+X3*(5.175290e-06)+X4*(-6.011153e-04)+
  X9*(2.637086e-03)+X10*(7.504337e-03)
detach(exp_clinical)
summary(exp_clinical$risk_score)
#divide the risk_score in to high or low risk by median
exp_clinical$score_ranking[exp_clinical$risk_score>=0.05919]=2
exp_clinical$score_ranking[exp_clinical$risk_score< 0.05919]=1
exp_clinical$score_ranking=factor(exp_clinical$score_ranking)
#head(exp_clinical)
#Kaplan meier analysis
fit.km=survfit(Surv(time,event)~score_ranking,data=exp_clinical)
#plot
ggsurvplot(fit.km, pval = TRUE, conf.int = F,
           risk.table = TRUE, risk.table.y.text.col = TRUE)
#summary
survdiff(Surv(time,event)~score_ranking,data=exp_clinical)


