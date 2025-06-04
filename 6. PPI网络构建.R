data<-read.table('string_node_degrees.tsv',quote='',header=F,row.names=1,check.names=F)
data$id<-rownames(data)
data<-data[,2:3]
names(data)[1]<-'degree'
data<-data[order(data$degree,decreasing = T),]
data<-data[data$degree>2,]
data<-data[order(data$degree),]
pdf('degree.pdf',width=5,height=6)
par(mar=c(2,6,3,2)+0.1)
barplot(data$degree,names.arg=data$id,col='blue',horiz =T,las=1,axis.lty=1,main='Gene degree')
dev.off()
