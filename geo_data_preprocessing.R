source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("affyPLM")
biocLite("RColorBrewer")
biocLite("limma")
biocLite("pheatmap")
###########################
#set the working directiory
setwd("E:\\geo_learn\\T")
setwd("E:/geo_learn/T")
library(affy)
Data<-ReadAffy()
sampleNames(Data)
image(Data[,1])
library(affyPLM)
library(RColorBrewer)
#convert the AffyBatch into PLMset
Pset<-fitPLM(Data)
#creat the corlor scheme
colors<-brewer.pal(5,"Set3")
#boxplot RLE(Relative log Expression)
Mbox(Pset,col=colors,main="RLE",las=1)
#assess RNA degradation in Affymetrix GeneChip data
data.deg<-AffyRNAdeg(Data)
#plot the RNA degradation
plotAffyRNAdeg(data.deg,col=colors)
#add legend
legend("topleft",sampleNames(Data),col = colors,lwd=1,inset = 0.05,cex = 0.5)
#Data cleaning
#converts the affybatch into expressionset rma(robust multi-array average)
eset.rma<-rma(Data)
#get the express data matrix from the expressionset
tumor_exprs<-exprs(eset.rma)
#preserved the dataset
write.csv(tumor_exprs,file = "tumor_exprs.csv")
###############################
#deal with the normal data
setwd("E:\\geo_learn\\N")
Data.normal<-ReadAffy()
sampleNames(Data.normal)
image(Data.normal[,1])
Pset<-fitPLM(Data.normal)
Mbox(Pset,col=colors,main="RLE",las=1)
data.deg<-AffyRNAdeg(Data.normal)
plotAffyRNAdeg(data.deg,col=colors)
legend("topleft",sampleNames(Data.normal),col = colors,lwd=1,inset = 0.05,cex = 0.5)
eset.rma<-rma(Data.normal)
normal_exprs<-exprs(eset.rma)
write.csv(normal_exprs,file = "normal_exprs.csv")
###############################
#combined the normal and cancer expression dataset
setwd("E:\\geo_learn\\combined")
#input the expression and information frome the GSE
expression.normal=read.csv(file = "normal_exprs.csv")
names(expression.normal)[1]="probe"
expression.tumor=read.csv(file = "tumor_exprs.csv")
names(expression.tumor)[1]="probe"
gene_names=read.table(file = "GPL97-17394_cleaned.txt",header = T,sep = '\t')
names(gene_names)[1]="probe"
#merge by the probe id
expression.whole=merge(expression.normal,expression.tumor,by="probe")
expression.whole=merge(expression.whole,gene_names,by="probe")
head(expression.whole)
