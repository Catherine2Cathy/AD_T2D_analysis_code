gene<-read.csv('hub_marker.csv',quote='',header=T,check.names = F)
tf<-read.table('trrust_rawdata.human.tsv',quote='',header=F,check.names=F)
tf<-tf[,1:2]
names(tf)<-c('id1','id2')
names(gene)<-'id1'
d1<-merge(tf,gene,by='id1')
names(gene)<-'id2'
d2<-merge(tf,gene,by='id2')
names(d2)<-c('id1','id2')
d<-rbind(d1,d2)
d<-d[!duplicated(d),]
names(d)<-c('node','target')
rna<-read.csv('RNA-RNA.csv',quote='',header=T,check.names=F)
rna<-read.csv('RNA-RNA.csv',quote='',header=T,check.names=F)
rna<-rna[,c(3,1,2)]
d$type<-'TF'
names(rna)[2]<-'target'
d<-rbind(d,rna)
write.csv(d,'interactopm.csv',quote=F)