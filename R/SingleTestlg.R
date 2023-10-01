pacman::p_load(expm,magic,glmnet,cluster,MASS,locfit,corpcor,caret,mvtnorm,pROC)
set.seed(seed1)
library(magic)
library(mvtnorm)
n_int<-224
p_X<-40
inv_r<-30
n_ext<-n_int*inv_r
#p_X<-10
p_A<-150
XAcolnames<-c(paste0("X",1:p_X),paste0("A",1:p_A))
intercept<- 0#-1/3 #-2.3
coef_X<-rep(0,p_X)
coef_A<-rep(0,p_A)
X_nonzero_index<-1:10
A_nonzero_index<-c(11,12,13,17,21,26,31,34,36,50,64,67,69,88,90)-10 #15 nonzero
coef_X[X_nonzero_index]<- 0.3/2
coef_A[A_nonzero_index]<- 0.25/2
coefXA<-c(coef_X,coef_A)
XAc<-list()
XAc[[10]]<-c(8,9)
XAc[[8]]<-c(1)
XAc[[6]]<-c(2)
XAc[[5]]<-c(3)
XAc[[4]]<-c(4)
XAc[[3]]<-c(10)
XAc[[2]]<-c(14)
XAc[[1]]<-c(11,13)

autocorM<-function(size,rho){return(rho^abs(row(diag(size))-col(diag(size))))}
generateM<-function(n,p_X,p_A,coefs,intercept,
                    rho=0.5,block_size=10,adj_r=0.3){
    Sigma_unit<-autocorM(block_size,rho)
    Sigma_list<-list()
    p<-p_X+p_A
    for(i in 1:(p/block_size)){Sigma_list[[i]]<-Sigma_unit}
    Sigma<-Reduce(adiag,Sigma_list)
    M<-mvtnorm::rmvnorm(n,sigma=Sigma)
    for(i in 1:length(XAc)){
        if(length(XAc[[i]])>0){
            for(j in XAc[[i]]){
                jj<-p_X+A_nonzero_index[j]
                M[,jj]<-M[,jj]+adj_r*M[,i]
            }
        }
    }
    M<-scale(M)
    y<-rbinom(n,size = 1,prob=expit(M%*%coefs+intercept))
    list("M"=M,"y"=y)
}
summary_stat<-function(M,y,p_X){
    study_info_scaled<-list()
    Xonlylm<-glm(y~.,data = data.frame(y,M[,1:p_X]),family = "binomial")
    study.m = list(Coeff=Xonlylm$coefficients[-1],
                   Covariance=vcov(Xonlylm)[-1,-1],Sample_size=nrow(M))
    study_info_scaled[[1]] <- study.m
    study_info_scaled_sep<-list()
    study_info_scaled_sep<-lapply(1:p_X, function(i){
        Xunilm<-glm(y~.,data = data.frame(y,M[,i]),family = "binomial")
        sum_Xunilm<-summary(Xunilm)
        study.m = list(Coeff=Xunilm$coefficients[2],
                       Covariance=vcov(Xunilm)[2,2],Sample_size=nrow(M),
                       P=sum_Xunilm$coefficients[2,4])
        study.m
    })
    return(list("multilm"=Xonlylm,"sum_mul"=study_info_scaled,"sum_uni"=study_info_scaled_sep))
}

My_ext<-generateM(n_ext,p_X,p_A,coefXA,intercept)
M_ext<-My_ext$M
y_ext<-My_ext$y
colnames(M_ext)<-XAcolnames
sum_res<-summary_stat(M_ext,y_ext,p_X)
multilm<-sum_res$multilm
sum_mul<-sum_res$sum_mul
sum_uni<-sum_res$sum_uni
# internal data
My_int<-generateM(n_int,p_X,p_A,coefXA,intercept)
M_int<-My_int$M
y_int<-c(My_int$y)
#y_int<-rbinom(nrow(M_int),size = 1,prob=expit(M_int%*%coefXA+intercept))
X_int<-M_int[,1:p_X]
A_int<-M_int[,-c(1:p_X)]
library(intgmm)
fuse_uni_1lam<-cv.intgmm(y_int,X_int,A_int,sum_uni,family = "binomial",penalty_type = "lasso",summary_type = "uni",tune_ratio = F,use_sparseC=T,shrink = T)
fuse_mul_1lam<-cv.intgmm(y_int,X_int,A_int,sum_mul,family = "binomial",penalty_type = "lasso",summary_type = "multi",tune_ratio = F,use_sparseC=T,shrink = T)
#fuse_mul_1lam1<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",tune_ratio = F,desparseC = T)

#M0<-readRDS("../FusionGMMdata/test2_40.rds")
#M0<-M0$M
#y0<-readRDS("../FusionGMMdata/test2binary_y_40.rds")
RC<-function(pred_y){
    suppressMessages(cur_auc<-c(pROC::auc(y0,expit(pred_y),direction = "<")))
}
fuse_lasso<-cv.glmnet(x = M_int,y = y_int,family = "binomial")
#suppressWarnings(r0<-RC(M0%*%coefXA))
#suppressWarnings(r0X<-RC(M0[,1:p_X]%*%coef_X))
#suppressWarnings(r0A<-RC(M0[,-c(1:p_X)]%*%coef_A))
#suppressWarnings(rX<-RC(M0[,1:p_X]%*%sum_mul[[1]]$Coeff))
#suppressWarnings(rlasso<-RC(c(M0%*%coef(fuse_lasso,s="lambda.min")[-1])))
#suppressWarnings(runi1<-RC(c(M0%*%fuse_uni_1lam$beta[-length(fuse_uni_1lam$beta)])))
#suppressWarnings(rmul1<-RC(c(M0%*%fuse_mul_1lam$beta[-length(fuse_uni_1lam$beta)])))
#a<-(round(c(rlasso,runi1,rmul1),3))
#a<-(round(c(r0,r0X,r0A,rX,rlasso,runi1,rmul1),3))
#names(a)<-c("rlasso","runi1","rmul1")
#print(a)
betamse<-c(sum((coef(fuse_lasso,s="lambda.min")[-1] - coefXA)^2),
sum((fuse_uni_1lam$beta[-length(fuse_uni_1lam$beta)] - coefXA)^2),
sum((fuse_mul_1lam$beta[-length(fuse_uni_1lam$beta)] - coefXA)^2))
names(betamse)<-c("rlasso","runi1","rmul1")
print(round(betamse,3))


fuse_uni_1lamF<-cv.intgmm(y_int,X_int,A_int,sum_uni,family = "binomial",penalty_type = "lasso",summary_type = "uni",tune_ratio = F,use_sparseC=T,shrink = F)
fuse_mul_1lamF<-cv.intgmm(y_int,X_int,A_int,sum_mul,family = "binomial",penalty_type = "lasso",summary_type = "multi",tune_ratio = F,use_sparseC=T,shrink = F)
betamse<-c(sum((coef(fuse_lasso,s="lambda.min")[-1] - coefXA)^2),
           sum((fuse_uni_1lamF$beta[-length(fuse_uni_1lamF$beta)] - coefXA)^2),
           sum((fuse_mul_1lamF$beta[-length(fuse_uni_1lamF$beta)] - coefXA)^2))
names(betamse)<-c("rlasso","runi1","rmul1")
print(round(betamse,3))
library(intgmm)
fuse_uni_1lam1<-cv.intgmm(scale(y_int,scale = F)/var(y_int),X_int,A_int,sum_uni,penalty_type = "lasso",summary_type = "uni",tune_ratio = F,use_sparseC=T)
fuse_mul_1lam1<-cv.intgmm(scale(y_int,scale = F)/var(y_int),X_int,A_int,sum_mul,penalty_type = "lasso",summary_type = "multi",tune_ratio = F,use_sparseC=T)
betamse<-c(sum((coef(fuse_lasso,s="lambda.min")[-1] - coefXA)^2),
           sum((fuse_uni_1lam1$beta - coefXA)^2),
           sum((fuse_mul_1lam1$beta- coefXA)^2))
names(betamse)<-c("rlasso","runi1","rmul1")
print(round(betamse,3))
print("end")

