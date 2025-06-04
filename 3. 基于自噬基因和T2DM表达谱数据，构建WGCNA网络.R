dataExpr<-read.csv('exprdata_all.csv',quote='',header=T,row.names=1,check.names=F)
gene<-read.csv('aut_gene.csv',quote='',header=T,check.names=F)
dataExpr<-dataExpr[rownames(dataExpr)%in%gene$id,]
traitData<-read.csv('phe_T2D.csv',quote='',header=T,row.names=1,check.names=F)
dataExpr<-as.data.frame(t(dataExpr))
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
pdf('sampletree.pdf',width=60,height=10)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
pdf('softthreshold.pdf',width=10,height=10)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()
power = sft$powerEstimate
k <- softConnectivity(datE=dataExpr,power=sft$powerEstimate)
k[which(is.na(k))]<-0
pdf('k.pdf',width=7,height=5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
TOMType = type, minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs=FALSE, corType = corType,
maxPOutliers=maxPOutliers, loadTOMs=TRUE,
verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
pdf('cluster.pdf',width=5,height=5)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
traitData$T2D<-ifelse(traitData$type=='T2D',1,0)
traitData$normal<-ifelse(traitData$type=='T2D',0,1)
traitData<-traitData[,-1:-2]
modTraitCor = cor(MEs_col, traitData, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf('correlation.pdf',width=6,height=6)
labeledHeatmap(Matrix = modTraitCor,
               xLabels = colnames(traitData),  
               yLabels = colnames(MEs_col),    
               cex.lab = 0.9,  
               ySymbols = colnames(MEs_col),  
               colorLabels = FALSE,   
               colors = blueWhiteRed(50),   
               textMatrix = textMatrix,   
               setStdMargins = TRUE,    
               cex.text = 1, zlim = c(-1,1),
               xLabelsAdj=1)  
dev.off()
x<-cbind(moduleLabels,moduleColors)
x<-as.data.frame(x)
data<-x[x$moduleColors=='blue'|x$moduleColors=='brown'|x$moduleColors=='yellow'|x$moduleColors=='green'|x$moduleColors=='red',]
data$name<-rownames(data)
data<-data[,-1]
dataExpr<-as.data.frame(t(dataExpr))
dataExpr$name<-rownames(dataExpr)
data<-merge(data,dataExpr,by='name')
rownames(data)<-data$name
data<-data[,-1:-2]
write.csv(data,'all_color.csv',quote=F)
