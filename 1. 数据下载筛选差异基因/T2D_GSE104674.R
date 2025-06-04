data<-read.table('GSE104674_rawCounts.txt.gz',quote='',header=T,row.names=1,check.names=F)
data$id<-gsub('\\:.*','',rownames(data))
data$id<-gsub('\\_.*','',data$id)
data<-aggregate(data[,-49],by=list(data$id),mean)
rownames(data)<-data$Group.1
data<-data[,-1]
library(GEOquery)
library(Biobase)
gse=getGEO("GSE104674",GSEMatrix = TRUE,destdir = ".",getGPL = T,AnnotGPL = T)
pdata<-pData(gse[[1]])
phe<-as.data.frame(names(data))
names(phe)<-'sample'
phe$type<-pdata$`subject status:ch1`
phe$type<-ifelse(phe$type=='Type 2 diabetes','T2D','normal')
library(edgeR)
y <- DGEList(counts=data,group=phe$type)
keep <- rowSums(cpm(y)>1) >= 2#至少在2个样本里cpm大于1
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
top<-topTags(et, n=nrow(data))
resdata<-as.data.frame(top)
deg<-subset(resdata, abs(resdata$logFC)>0.3 & FDR < 0.05)
resdata$group<-ifelse(resdata$FDR<0.05 & abs(resdata$logFC)>0.3,ifelse(resdata$logFC>0.3,'up','down'),'no')
library(ggplot2)
p<-ggplot(resdata,aes(x=logFC,y=-log10(FDR)))
p<-p+geom_point(aes(col=group))+scale_color_manual(values=c(up='red',down='blue',no='grey'))
p<-p+geom_hline(yintercept = -log10(0.05), lty=4,lwd=0.6)+geom_vline(xintercept = c(-0.3, 0.3), lty=4,lwd=0.6)
p <- p+labs(x=expression(log2FoldChange),y=expression(-log10(FDR)))+theme_bw()+xlim(-7,7)
ggsave('volcanoplot.pdf',width=5,height=5)
newdata<-as.data.frame(y$pseudo.counts)
newdata$id<-rownames(newdata)
newdata<-newdata[,-49]
newdata<-log2(newdata+1)
write.csv(newdata,'exprdata_all.csv',quote=F)
write.csv(phe,'phe_T2D.csv',quote=F)
write.csv(deg,'deg_T2D.csv',quote=F)





