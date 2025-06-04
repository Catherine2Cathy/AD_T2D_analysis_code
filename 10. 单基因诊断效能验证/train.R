d1<-read.csv('markers_AD.csv',quote='',header=T,row.names=1,check.names=F)
d2<-read.csv('markers_T2DM.csv',quote='',header=T,row.names=1,check.names=F)
p1<-read.csv('phe_ad.csv',quote='',header=T,row.names=1,check.names=F)
p2<-read.csv('phe_T2D.csv',quote='',header=T,row.names=1,check.names=F)
hub<-read.csv('hub_marker.csv',quote='',header=T)
d1<-d1[rownames(d1)%in%hub$id,]
d2<-d2[rownames(d2)%in%hub$id,]
d1<-as.data.frame(t(d1))
d2<-as.data.frame(t(d2))
d1<-cbind(p1,d1)
d2<-cbind(p2,d2)
d1$type_reverse<-ifelse(d1$type=='AD',0,1)
d2$type_reverse<-ifelse(d2$type=='T2D',0,1)
d1$type<-ifelse(d1$type=='AD',1,0)
d2$type<-ifelse(d2$type=='T2D',1,0)

pacman::p_load(tidyverse,pROC,verification,ROCR)
run_home <- "./"
gene_list <- c("ANXA5","BAG3","CDKN1A","MET","PFKFB3")

d1$type <- as.factor(d1$type)
p_list <- map(gene_list[-4],function(gene){
    pred <- prediction(d1[,gene], d1$type)
    auc_min <- performance(pred,"auc")@y.values[[1]]
    perf_min <- performance(pred,"tpr","fpr")
    roc_res <- roc(type ~ get0(gene),data=d1)
    ci_info <- ci.auc(roc_res)
    labels <- as.numeric(d1$type) - 1
    p_value <- roc.area(labels,roc_res$predictor)$p.value
    auc_ci <- paste0(gene,"_AUC = ",round(auc_min,3),"\n",sprintf("95%% CI: %.3f-%.3f", ci_info[1], ci_info[3]),"\n",sprintf("p_value: %0.2e", p_value))
    print(auc_ci)
    pdf(str_glue("{run_home}/results/train_ad_{gene}_roc.pdf"),width=5,height=5.5)
    plot(perf_min,colorize=FALSE, col="red")
    text(0.8,0.2, labels = auc_ci)
    lines(c(0,1),c(0,1), lty = 4)
    dev.off()
})

d2$type <- as.factor(d2$type)
p_list <- map(gene_list[-c(4,5)],function(gene){
    pred <- prediction(d2[,gene], d2$type)
    auc_min <- performance(pred,"auc")@y.values[[1]]
    perf_min <- performance(pred,"tpr","fpr")
    roc_res <- roc(type ~ get0(gene),data=d2)
    ci_info <- ci.auc(roc_res)
    labels <- as.numeric(d2$type) - 1
    p_value <- roc.area(labels,roc_res$predictor)$p.value
    auc_ci <- paste0(gene,"_AUC = ",round(auc_min,3),"\n",sprintf("95%% CI: %.3f-%.3f", ci_info[1], ci_info[3]),"\n",sprintf("p_value: %0.2e", p_value))
    print(auc_ci)
    pdf(str_glue("{run_home}/results/train_t2d_{gene}_roc.pdf"),width=5,height=5.5)
    plot(perf_min,colorize=FALSE, col="red")
    text(0.8,0.2, labels = auc_ci)
    lines(c(0,1),c(0,1), lty = 4)
    dev.off()
})

d1$type_reverse <- as.factor(d1$type_reverse)
p_list <- map(gene_list[4],function(gene){
    pred <- prediction(d1[,gene], d1$type_reverse)
    auc_min <- performance(pred,"auc")@y.values[[1]]
    perf_min <- performance(pred,"tpr","fpr")
    roc_res <- roc(type_reverse ~ get0(gene),data=d1)
    ci_info <- ci.auc(roc_res)
    labels <- as.numeric(d1$type_reverse) - 1
    p_value <- roc.area(labels,roc_res$predictor)$p.value
    auc_ci <- paste0(gene,"_AUC = ",round(auc_min,3),"\n",sprintf("95%% CI: %.3f-%.3f", ci_info[1], ci_info[3]),"\n",sprintf("p_value: %0.2e", p_value))
    print(auc_ci)
    pdf(str_glue("{run_home}/results/train_ad_{gene}_roc.pdf"),width=5,height=5.5)
    plot(perf_min,colorize=FALSE, col="red")
    text(0.8,0.2, labels = auc_ci)
    lines(c(0,1),c(0,1), lty = 4)
    dev.off()
})

d2$type_reverse <- as.factor(d2$type_reverse)
p_list <- map(gene_list[c(4,5)],function(gene){
    pred <- prediction(d2[,gene], d2$type_reverse)
    auc_min <- performance(pred,"auc")@y.values[[1]]
    perf_min <- performance(pred,"tpr","fpr")
    roc_res <- roc(type_reverse ~ get0(gene),data=d2)
    ci_info <- ci.auc(roc_res)
    labels <- as.numeric(d2$type_reverse) - 1
    p_value <- roc.area(labels,roc_res$predictor)$p.value
    auc_ci <- paste0(gene,"_AUC = ",round(auc_min,3),"\n",sprintf("95%% CI: %.3f-%.3f", ci_info[1], ci_info[3]),"\n",sprintf("p_value: %0.2e", p_value))
    print(auc_ci)
    pdf(str_glue("{run_home}/results/train_t2d_{gene}_roc.pdf"),width=5,height=5.5)
    plot(perf_min,colorize=FALSE, col="red")
    text(0.8,0.2, labels = auc_ci)
    lines(c(0,1),c(0,1), lty = 4)
    dev.off()
})
