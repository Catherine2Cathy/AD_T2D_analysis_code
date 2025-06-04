da1<-read.csv('deg_AD.csv',quote='',header=T)
da2<-read.csv('deg_T2D.csv',quote='',header=T)
da3<-read.csv('turquoise+yellow+green.csv',quote='',header=T)
da4<-read.csv('all_color.csv',quote='',header=T)
d1<-da1[,1]
d2<-da2[,1]
d3<-da3[,1]
d4<-da4[,1]
library(VennDiagram)
v<-venn.diagram(x=list(DEG_AD=d1,DEG_T2DM=d2,WGCNA_AD=d3,WGCNA_T2DM=d4),col="transparent",fill=c("red", "blue", "green","gold"),filename=NULL,alpha =0.4,cat.cex = 1.5,cat.fontfamily = "serif",cat.default.pos = "text",cat.col = c("darkred", "darkblue", "darkgreen","gold"))
pdf(file="all.pdf")
grid.draw(v)
dev.off()

d1<-as.data.frame(d1)
d2<-as.data.frame(d2)
d3<-as.data.frame(d3)
d4<-as.data.frame(d4)
names(d1)<-'id'
names(d2)<-'id'
names(d3)<-'id'
names(d4)<-'id'
d<-merge(merge(d1,d2,by='id'),merge(d3,d4,by='id'),by='id')
write.csv(d,'marker.csv',quote=F)

da4<-da4[da4$X%in%d$id,]
da3<-da3[da3$X%in%d$id,]
rownames(da3)<-da3$X
da3<-da3[,-1]
rownames(da4)<-da4$X
da4<-da4[,-1]
write.csv(da3,'markers_AD.csv',quote=F)
write.csv(da4,'markers_T2DM.csv',quote=F)