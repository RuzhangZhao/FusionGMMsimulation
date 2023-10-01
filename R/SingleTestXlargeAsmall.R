
n_int<-500
inv_r<-10
p_X<-10
print(paste0("p_X=",p_X,";inv_r=",inv_r,":n_int=",n_int))
#source("/users/rzhao1/fusion/code/fusionGMM.R")
library(magic)
library(mvtnorm)
#n_int<-100
#inv_r<-3
n_ext<-n_int*inv_r
#p_X<-10
p_A<-150
intercept<- -1.6
coef_X<-rep(0,p_X)
coef_A<-rep(0,p_A)
X_nonzero_index<-1:10
A_nonzero_index<-c(11,12,13,17,21,26,31,34,36,50,64,67,69,88,90)-10 #15 nonzero
coef_X[X_nonzero_index]<- 0.3/3
coef_A[A_nonzero_index]<- 0.25/3
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
    y<-M%*%coefs+intercept+rnorm(n,0,3/3)
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

R2<-function(y0,predy){
    1-sum((predy-scale(y0,scale = F))^2)/sum((scale(y0,scale = F))^2)
}
rr_gap<-list()

My_ext<-generateM(n_ext,p_X,p_A,coefXA,intercept)
M_ext<-My_ext$M
y_ext<-My_ext$y
M_ext<-M_ext[,c((p_X+1):(p_X+p_A),1:p_X)]
AXcolnames<-c(paste0("X",1:p_A),paste0("A",1:p_X))

colnames(M_ext)<-AXcolnames
sum_res<-summary_stat(M_ext,y_ext,p_A)
multilm<-sum_res$multilm
sum_mul<-sum_res$sum_mul
sum_uni<-sum_res$sum_uni
# internal data
My_int<-generateM(n_int,p_X,p_A,coefXA,intercept)
M_int<-My_int$M
M_int<-M_int[,c((p_X+1):(p_X+p_A),1:p_X)]
y_int<-c(My_int$y)
y_int<-scale(y_int,scale = F)
X_int<-M_int[,1:p_A]
A_int<-M_int[,-c(1:p_A)]

fuse_uni_1lam<-fusionGMM.addition(X_int,A_int,y_int,sum_uni,penalty_type = "lasso",summary_type = "uni",approx_cross_validation =F,tune_ratio = F)
fuse_mul_1lam<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = F)

fuse_lasso<-cv.glmnet(x = M_int,y = y_int)

r0<-R2(y_ext,M_ext%*%coefAX)
rlasso<-R2(y_ext,M_ext%*%coef(fuse_lasso,s="lambda.min")[-1])

runi1<-R2(y_ext,M_ext%*%fuse_uni_1lam$beta)
rmul1<-R2(y_ext,M_ext%*%fuse_mul_1lam$beta)
print(c(r0,rlasso,runi1,rmul1))
