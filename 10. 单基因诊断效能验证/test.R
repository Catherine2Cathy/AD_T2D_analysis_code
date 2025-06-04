d1<-read.csv('AD.csv',quote='',header=T,row.names=1,check.names=F)
d2<-read.csv('T2D.csv',quote='',header=T,row.names=1,check.names=F)
d1$type_reverse<-ifelse(d1$type=='AD',0,1)
d2$type_reverse<-ifelse(d2$type=='T2D',0,1)
d1$type<-ifelse(d1$type=='AD',1,0)
d2$type<-ifelse(d2$type=='T2D',1,0)

library(ROCR)
re11<-d1[,c(2,1)]
names(re11)<-c('fp','event')
re11$event=as.factor(re11$event)
pred11<- prediction(re11[,1], re11[,2])
auc_min11 = performance(pred11,"auc")@y.values[[1]]
perf_min11 <- performance(pred11,"tpr","fpr")
pdf('roc11.pdf',width=5,height=5.5)
plot(perf_min11,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("ANXA5_AUC = ",round(auc_min11,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re12<-d1[,c(3,1)]
names(re12)<-c('fp','event')
re12$event=as.factor(re12$event)
pred12<- prediction(re12[,1], re12[,2])
auc_min12 = performance(pred12,"auc")@y.values[[1]]
perf_min12 <- performance(pred12,"tpr","fpr")
pdf('roc12.pdf',width=5,height=5.5)
plot(perf_min12,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("BAG3_AUC = ",round(auc_min12,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re13<-d1[,c(7,1)]
names(re13)<-c('fp','event')
re13$event=as.factor(re13$event)
pred13<- prediction(re13[,1], re13[,2])
auc_min13 = performance(pred13,"auc")@y.values[[1]]
perf_min13 <- performance(pred13,"tpr","fpr")
pdf('roc13.pdf',width=5,height=5.5)
plot(perf_min13,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("CDKN1A_AUC = ",round(auc_min13,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re14<-d1[,c(11,14)]
names(re14)<-c('fp','event')
re14$event=as.factor(re14$event)
pred14<- prediction(re14[,1], re14[,2])
auc_min14 = performance(pred14,"auc")@y.values[[1]]
perf_min14 <- performance(pred14,"tpr","fpr")
pdf('roc14.pdf',width=5,height=5.5)
plot(perf_min14,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("MET_AUC = ",round(auc_min14,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re15<-d1[,c(12,1)]
names(re15)<-c('fp','event')
re15$event=as.factor(re15$event)
pred15<- prediction(re15[,1], re15[,2])
auc_min15 = performance(pred15,"auc")@y.values[[1]]
perf_min15 <- performance(pred15,"tpr","fpr")
pdf('roc15.pdf',width=5,height=5.5)
plot(perf_min15,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("PFKFB3_AUC = ",round(auc_min15,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re21<-d2[,c(2,1)]
names(re21)<-c('fp','event')
re21$event=as.factor(re21$event)
pred21<- prediction(re21[,1], re21[,2])
auc_min21 = performance(pred21,"auc")@y.values[[1]]
perf_min21 <- performance(pred21,"tpr","fpr")
pdf('roc21.pdf',width=5,height=5.5)
plot(perf_min21,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("ANXA5_AUC = ",round(auc_min21,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re22<-d2[,c(3,1)]
names(re22)<-c('fp','event')
re22$event=as.factor(re22$event)
pred22<- prediction(re22[,1], re22[,2])
auc_min22 = performance(pred22,"auc")@y.values[[1]]
perf_min22 <- performance(pred22,"tpr","fpr")
pdf('roc22.pdf',width=5,height=5.5)
plot(perf_min22,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("BAG3_AUC = ",round(auc_min22,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re23<-d2[,c(7,1)]
names(re23)<-c('fp','event')
re23$event=as.factor(re23$event)
pred23<- prediction(re23[,1], re23[,2])
auc_min23 = performance(pred23,"auc")@y.values[[1]]
perf_min23 <- performance(pred23,"tpr","fpr")
pdf('roc23.pdf',width=5,height=5.5)
plot(perf_min23,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("CDKN1A_AUC = ",round(auc_min23,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re24<-d2[,c(10,13)]
names(re24)<-c('fp','event')
re24$event=as.factor(re24$event)
pred24<- prediction(re24[,1], re24[,2])
auc_min24 = performance(pred24,"auc")@y.values[[1]]
perf_min24 <- performance(pred24,"tpr","fpr")
pdf('roc24.pdf',width=5,height=5.5)
plot(perf_min24,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("MET_AUC = ",round(auc_min24,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()

re25<-d2[,c(11,13)]
names(re25)<-c('fp','event')
re25$event=as.factor(re25$event)
pred25<- prediction(re25[,1], re25[,2])
auc_min25 = performance(pred25,"auc")@y.values[[1]]
perf_min25 <- performance(pred25,"tpr","fpr")
pdf('roc25.pdf',width=5,height=5.5)
plot(perf_min25,colorize=FALSE, col="red")
text(0.8,0.2, labels = paste0("PFKFB3_AUC = ",round(auc_min25,3)))
lines(c(0,1),c(0,1), lty = 4 )
dev.off()