d1<-read.csv('markers_AD.csv',quote='',header=T,row.names=1,check.names=F)
d2<-read.csv('markers_T2DM.csv',quote='',header=T,row.names=1,check.names=F)
hub<-read.csv('hub_marker.csv',quote='',header=T)
d1<-d1[rownames(d1)%in%hub$id,]
d2<-d2[rownames(d2)%in%hub$id,]
p1<-read.csv('phe_ad.csv',quote='',header=T,row.names=1,check.names=F)
p2<-read.csv('phe_T2D.csv',quote='',header=T,row.names=1,check.names=F)
d1<-as.data.frame(t(d1))
d2<-as.data.frame(t(d2))
d1<-cbind(p1,d1)
d2<-cbind(p2,d2)

d1<-d1[,-1]
d2<-d2[,-1]
library(reshape2)
c1<-melt(d1,id.vars='type')
library(ggsignif)
library(ggforce)
library(ggpubr)
pdf('AD.pdf',width=6,height=4)
ggboxplot(c1, x="variable", y="value",outlier.shape = NA,
          fill = "type",xlab = "", ylab = "expression level") +
    stat_compare_means(aes(group= type),label = "p.signif",label.y = 14,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", " "))) +
    scale_fill_manual(values =c('red','blue')) +
    theme_classic()+
    #theme(legend.position = "none")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
c2<-melt(d2,id.vars='type')
c2$type<-factor(c2$type,levels = c('T2D','normal'))
pdf('T2DM.pdf',width=6,height=4)
ggboxplot(c2, x="variable", y="value",outlier.shape = NA,
          fill = "type",xlab = "", ylab = "expression level") +
    stat_compare_means(aes(group= type),label = "p.signif",label.y = 13,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", " "))) +
    scale_fill_manual(values =c('red','blue')) +
    theme_classic()+
    #theme(legend.position = "none")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()