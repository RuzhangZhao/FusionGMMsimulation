source("/users/rzhao1/intgmm/code/intgmm.R")
#path<-"/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/nonpos/"
path<-"/users/rzhao1/intgmm/result/nonpos/"
library(magic)
library(mvtnorm)
p_X<-40
inv_r<-39
p_A<-150
XAcolnames<-c(paste0("X",1:p_X),paste0("A",1:p_A))
intercept<- -1/3
coef_X<-rep(0,p_X)
coef_A<-rep(0,p_A)
X_nonzero_index<-1:10
A_nonzero_index<-c(11,12,13,17,21,26,31,34,36,50,64,67,69,88,90)-10 #15 nonzero
coef_X[X_nonzero_index]<- 0.3/2
coef_A[A_nonzero_index]<- 0.25/2
coefXA<-c(coef_X,coef_A)
coefAX<-c(coef_A,coef_X)
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
    M<-rmvnorm(n,sigma=Sigma)
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

### Infinity Test Data for Validation
# first 150 vs last 40
if(F){
    ntest = 40000
    My0<-generateM(ntest,p_X,p_A,coefXA,intercept)
    M0<-My0$M
    saveRDS(M0,paste0("/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/testonlyX_",p_X,".rds"))
}else{
    #M0<-readRDS(paste0("/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/testonlyX_",p_X,".rds"))
    M0<-readRDS(paste0("/users/rzhao1/intgmm/data/fdrdata/testonlyX_",p_X,".rds"))
    M0<-M0[,c((1+p_X):(p_A+p_X),1:p_X)]
}
library(pROC)
RC<-function(pred_y){
    suppressMessages(cur_auc<-c(auc(y0,expit(pred_y),direction = "<")))
}

internal_index<-c(1:floor(nrow(M0)/(1+inv_r)))
M_int<-M0[internal_index,]
M_ext<-M0[-internal_index,]
correct_pos<-which(coefAX!=0)
nonpostotal<-list()
for(i in 1:100){
    message(paste0("EPOCH",i))
    y0<-rbinom(nrow(M0),size = 1,prob=expit(M0%*%coefAX+intercept))
    print(RC(M0%*%coefAX))
    y_int<-y0[internal_index]
    y_ext<-y0[-internal_index]
    sum_res<-summary_stat(M_ext,y_ext,p_A)
    multilm<-sum_res$multilm
    sum_mul<-sum_res$sum_mul
    sum_uni<-sum_res$sum_uni

    X_int<-M_int[,1:p_A]
    A_int<-M_int[,-c(1:p_A)]

    fuse_uni_ada<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",inference = T)
    fuse_uni_ada_1lam<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",tune_ratio = F,inference = T)

    fuse_mul_ada<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",inference = T)
    fuse_mul_ada_1lam<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",tune_ratio = F,inference = T)

    if(F){
        #fuse_lasso<-cv.glmnet(x = M_int,y = y_int)
        fuse_ridge<-cv.glmnet(x = M_int,y = y_int,alpha=0)
        www<-1/abs(c(coef(fuse_ridge,s="lambda.min")[-1]))^(1/2)
        www[is.infinite(www)]<-max(www[!is.infinite(www)])*100
        fuse_ada<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www)
    }
    fdrgmm<-function(fuse_uni_ada){
        sum(!fuse_uni_ada[['corrected_pos']]%in%correct_pos)/length(fuse_uni_ada[['corrected_pos']])
    }
    power<-function(fuse_uni_ada){
        sum(fuse_uni_ada[['corrected_pos']]%in%correct_pos)/length(correct_pos)
    }
    nonpos<-c(fdrgmm(fuse_uni_ada),fdrgmm(fuse_uni_ada_1lam),fdrgmm(fuse_mul_ada),fdrgmm(fuse_mul_ada_1lam))
    nonpos<-c(nonpos,c(power(fuse_uni_ada),power(fuse_uni_ada_1lam),power(fuse_mul_ada),power(fuse_mul_ada_1lam)))
    names(nonpos)<-c("fdr_uni_ada","fdr_uni_ada_1lam","fdr_mul_ada","fdr_mul_ada_1lam","power_uni_ada","power_uni_ada_1lam","power_mul_ada","power_mul_ada_1lam")
    nonpostotal[[i]]<-nonpos
    print(nonpos)
    if (i %% 5 == 0){
        saveRDS(nonpostotal,paste0(path,'binary_pos_ir_',inv_r,'_pX_',p_X,'.rds'))
    }
}
