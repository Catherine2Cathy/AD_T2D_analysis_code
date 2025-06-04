library(GEOquery)
library(Biobase)
gse=getGEO('GSE122063',GSEMatrix = TRUE,destdir = ".",getGPL = T,AnnotGPL = T)
exprs<-exprs(gse[[1]])
pdata<-pData(gse[[1]])
fdata<-fData(gse[[1]])
data<-as.data.frame(exprs)
data$ID<-rownames(data)
gene_id<-fdata[,c(1,10)]
names(gene_id)<-c('ID','symbol')
data<-merge(gene_id,data,by='ID')
data<-data[,-1]
gene<-read.csv('hub_marker.csv',quote='',header=T)
data<-data[data$symbol%in%gene$id,]
data<-aggregate(data[,-1],by=list(data$symbol),mean)
rownames(data)<-data$Group.1
data<-data[,-1]
data<-log2(data+1)
data<-as.data.frame(t(data))
phe<-pdata[c(2,43)]
names(phe)<-c('sample','type')
phe<-phe[phe$type!='Vascular dementia',]
phe$type<-ifelse(phe$type=="Alzheimer's disease",'AD','normal')
data$sample<-rownames(data)
data<-merge(phe,data,by='sample')
rownames(data)<-data$sample
data<-data[,-1]
write.csv(data,'AD.csv',quote=F)

library(reshape2)
c1<-melt(data,id.vars='type')
library(ggsignif)
library(ggforce)
library(ggpubr)
pdf('AD.pdf',width=6,height=4)
ggboxplot(c1, x="variable", y="value",outlier.shape = NA,
          fill = "type",xlab = "", ylab = "expression level") +
    stat_compare_means(aes(group= type),label = "p.signif",label.y = 5,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) +
    scale_fill_manual(values =c('red','blue')) +
    theme_classic()+
    #theme(legend.position = "none")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
