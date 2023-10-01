a<-readRDS("../FusionGMMdata/nonpos/pos_ir_39_pX_40.rds")
truepos<-c(3,7,26,151,152,159)
t1e<-rep(0,5)
fdr<-rep(0,5)
for(j in 1:100){
    for(i in 1:5){
        curpos<-a[[j]][[i]]
        t1e[i]<-t1e[i]+sum(!curpos%in%truepos)/c(sum(!curpos%in%truepos)+190-length(union(curpos,truepos)))
        fdr[i]<-fdr[i]+sum(!curpos%in%truepos)/length(curpos)
    }
}
names(t1e)<-c("adalasso","intGMM-ada-uni-2lambda","intGMM-ada-uni-1lambda","intGMM-ada-mul-2lambda","intGMM-ada-mul-1lambda")
names(fdr)<-c("adalasso","intGMM-ada-uni-2lambda","intGMM-ada-uni-1lambda","intGMM-ada-mul-2lambda","intGMM-ada-mul-1lambda")

fdr<-c()
pwr<-c()
fdr2<-c()
pwr2<-c()
coverage<-c()
for(i in 1:length(res)){
    fdr<-rbind(fdr,res[[i]][[1]])
    fdr2<-rbind(fdr2,res[[i]][[2]])
    pwr<-rbind(pwr,res[[i]][[3]])
    pwr2<-rbind(pwr2,res[[i]][[4]])
    coverage<-rbind(coverage,res[[i]][[5]])
}
colMeans(fdr)
colMeans(pwr)
colMeans(coverage)

colMeans(fdr2)
colMeans(pwr2)

