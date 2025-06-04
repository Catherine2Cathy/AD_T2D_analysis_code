library(GEOquery)
library(Biobase)
gse=getGEO("GSE109887",GSEMatrix = TRUE,destdir = ".",getGPL = T,AnnotGPL = T)
exprs<-exprs(gse[[1]])
pdata<-pData(gse[[1]])
data<-as.data.frame(exprs)
phe<-pdata[c(2,34)]
names(phe)<-c('sample','type')
phe$type<-ifelse(phe$type=='AD','AD','normal')
data<-as.data.frame(t(data))
data$sample<-rownames(data)
data<-merge(phe,data,by='sample')
rownames(data)<-data$sample
phe<-data[,1:2]
data<-data[,-1:-2]
data<-as.data.frame(t(data))
write.csv(data,'exprdata.csv',quote=F)
write.csv(phe,'phe_ad.csv',quote=F)

library(limma)
group<- factor(phe$type,levels = c('AD','normal'))
design <- model.matrix(~0+group)
colnames(design) <- c('AD','normal')
rownames(design) <- rownames(phe)
fit <- lmFit(data,design)
contrast.matrix <- makeContrasts(AD - normal,levels=design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2<- eBayes(fit2)
all_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 30000,sort.by = 'logFC')
deg<- subset(all_diff,abs(all_diff$logFC)>0.3 &all_diff$adj.P.Val<0.05)
resdata<-all_diff[,2:6]
resdata$group<-ifelse(resdata$adj.P.Val<0.05 & abs(resdata$logFC)>0.3,ifelse(resdata$logFC>0.3,'up','down'),'no')
library(ggplot2)
p<-ggplot(resdata,aes(x=logFC,y=-log10(adj.P.Val)))
p<-p+geom_point(aes(col=group))+scale_color_manual(values=c(up='red',down='blue',no='grey'))
p<-p+geom_hline(yintercept = -log10(0.05), lty=4,lwd=0.6)+geom_vline(xintercept = c(-0.3, 0.3), lty=4,lwd=0.6)
p<-p+theme_bw()
pdf('volcanoplot.pdf',width=5,height=5)
p
dev.off()
write.csv(deg,'deg_AD.csv',quote=F)

