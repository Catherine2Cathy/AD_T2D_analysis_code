library(GEOquery)
library(Biobase)
gse=getGEO('GSE64998',GSEMatrix = TRUE,destdir = ".",getGPL = T,AnnotGPL = T)
exprs<-exprs(gse[[1]])
pdata<-pData(gse[[1]])
fdata<-fData(gse[[1]])
data<-as.data.frame(exprs)
data$ID<-rownames(data)
gene_id<-fdata[,c(1,3)]
names(gene_id)<-c('ID','symbol')
data<-merge(gene_id,data,by='ID')
data<-data[,-1]
gene<-read.csv('hub_marker.csv',quote='',header=T)
data<-data[data$symbol%in%gene$id,]
data<-aggregate(data[,-1],by=list(data$symbol),max)
rownames(data)<-data$Group.1
data<-data[,-1]
data<-as.data.frame(t(data))
phe<-pdata[c(2,34)]
names(phe)<-c('sample','type')
data$sample<-rownames(data)
data<-merge(phe,data,by='sample')
rownames(data)<-data$sample
data<-data[,-1]
data$type<-ifelse(data$type=='obese type 2 diabetic patient','T2D','normal')
write.csv(data,'T2D.csv',quote=F)


library(reshape2)
c1<-melt(data,id.vars='type')
library(ggsignif)
library(ggforce)
library(ggpubr)
c1$type<-factor(c1$type,levels = c('T2D','normal'))
pdf('T2D.pdf',width=6,height=4)
ggboxplot(c1, x="variable", y="value",outlier.shape = NA,
          fill = "type",xlab = "", ylab = "expression level") +
    stat_compare_means(aes(group= type),label = "p.signif",label.y = 12,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))) +
    scale_fill_manual(values =c('red','blue')) +
    theme_classic()+
    #theme(legend.position = "none")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
