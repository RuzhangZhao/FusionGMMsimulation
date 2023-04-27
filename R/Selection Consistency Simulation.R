library(magic)
library(mvtnorm)
p_X<-40
inv_r<-39
p_A<-150
XAcolnames<-c(paste0("X",1:p_X),paste0("A",1:p_A))
intercept<- -2.3
coef_X<-rep(0,p_X)
coef_A<-rep(0,p_A)
X_nonzero_index<-1:10
A_nonzero_index<-c(11,12,13,17,21,26,31,34,36,50,64,67,69,88,90)-10 #15 nonzero
coef_X[X_nonzero_index]<- 0.3
coef_A[A_nonzero_index]<- 0.25
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
    M<-rmvnorm(n,sigma=Sigma)
    for(i in 1:length(XAc)){
        if(length(XAc[[i]])>0){
            for(j in XAc[[i]]){
                jj<-p_X+A_nonzero_index[j]
                M[,jj]<-M[,jj]+adj_r*M[,i]
            }
        }
    }
    M<-M[,c((p_X+1):(p_X+p_A),1:p_X)]
    M<-scale(M)
    y<-M%*%coefs+intercept+rnorm(n,0,3)
    list("M"=M,"y"=y)
}
summary_stat<-function(M,y,p_X){
    study_info_scaled<-list()
    Xonlylm<-lm(y~.,data = data.frame(y,M[,1:p_X]))
    study.m = list(Coeff=Xonlylm$coefficients[-1],
                   Covariance=vcov(Xonlylm)[-1,-1],Sample_size=nrow(M))
    study_info_scaled[[1]] <- study.m
    study_info_scaled_sep<-list()
    study_info_scaled_sep<-lapply(1:p_X, function(i){
        Xunilm<-lm(y~.,data = data.frame(y,M[,i]))
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
    ntest = 10^5
    My0<-generateM(ntest,p_X,p_A,coefXA,intercept)
    M0<-My0$M
    saveRDS(M0,paste0("/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/testonlyX_",p_X,".rds"))
}else{
    M0<-readRDS(paste0("/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/testonlyX_",p_X,".rds"))
}

R2<-function(y0,predy){
    1-sum((predy-scale(y0,scale = F))^2)/sum((scale(y0,scale = F))^2)
}

X_select<-c(1,2,9)
A_select<-c(3,4,9)
nonzero_index<-c(A_nonzero_index[A_select],p_A+X_nonzero_index[X_select])
coefXA<-rep(0,p_X+p_A)
nonpostotal<-list()
for(i in 1:100){
    coefXA[nonzero_index]<-rnorm(length(nonzero_index),0.1,0.01)
    y0<-c(M0%*%coefXA+rnorm(nrow(M0),0,1))
    R2(y0,M0%*%coefXA)
    internal_index<-c(1:floor(nrow(M0)/(1+inv_r)))
    M_int<-M0[internal_index,]
    y_int<-y0[internal_index]
    M_ext<-M0[-internal_index,]
    y_ext<-y0[-internal_index]
    sum_res<-summary_stat(M_ext,y_ext,p_A)
    multilm<-sum_res$multilm
    sum_mul<-sum_res$sum_mul
    sum_uni<-sum_res$sum_uni

    X_int<-M_int[,1:p_A]
    A_int<-M_int[,-c(1:p_A)]

    fuse_uni_ada<-fusionGMM.addition(X_int,A_int,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",approx_cross_validation =F,inference = T)
    fuse_uni_ada_1lam<-fusionGMM.addition(X_int,A_int,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",approx_cross_validation =F,tune_ratio = F,inference = T)

    fuse_mul_ada<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",approx_cross_validation =F,inference = T)
    fuse_mul_ada_1lam<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = F,inference = T)

    #fuse_lasso<-cv.glmnet(x = M_int,y = y_int)
    fuse_ridge<-cv.glmnet(x = M_int,y = y_int,alpha=0)
    www<-1/abs(c(coef(fuse_ridge,s="lambda.min")[-1]))^(1/2)
    www[is.infinite(www)]<-max(www[!is.infinite(www)])*100
    fuse_ada<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www)
    nonpos<-list()
    nonpos[[1]]<-which(coef(fuse_ada,s = "lambda.1se")[-1]!=0)
    nonpos[[2]]<-fuse_uni_ada$corrected_pos
    nonpos[[3]]<-fuse_uni_ada_1lam$corrected_pos
    nonpos[[4]]<-fuse_mul_ada$corrected_pos
    nonpos[[5]]<-fuse_mul_ada_1lam$corrected_pos
    nonpostotal[[i]]<-nonpos
    if (i %% 5 == 0){
        saveRDS(nonpostotal,paste0('/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/nonpos/pos_ir_',inv_r,'_pX_',p_X,'.rds'))
    }
}


