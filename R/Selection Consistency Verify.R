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
t1e<-t1e/100
fdr<-fdr/100
names(t1e)<-c("adalasso","intGMM-ada-uni-2lambda","intGMM-ada-uni-1lambda","intGMM-ada-mul-2lambda","intGMM-ada-mul-1lambda")
names(fdr)<-c("adalasso","intGMM-ada-uni-2lambda","intGMM-ada-uni-1lambda","intGMM-ada-mul-2lambda","intGMM-ada-mul-1lambda")

> t1e
adalasso intGMM-ada-uni-2lambda intGMM-ada-uni-1lambda
1.141304e-03           2.065217e-03           2.500000e-03
intGMM-ada-mul-2lambda intGMM-ada-mul-1lambda
5.434783e-05           5.434783e-05
> fdr
adalasso intGMM-ada-uni-2lambda intGMM-ada-uni-1lambda
0.030452381            0.055309524            0.066666667
intGMM-ada-mul-2lambda intGMM-ada-mul-1lambda
0.001428571            0.001428571
