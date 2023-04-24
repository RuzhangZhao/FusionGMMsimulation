

id1<-1
id2<-2
id3<-3
p_X_list<-c(10,40)
#inv_r_list<-c(3,6,10,30)
inv_r_list<-c(30,10,6,3)
n_int_list<-floor(exp(seq(log(100),log(3000),(log(3000) - log(100))/9)))
n_int_list<-floor((seq((100),(3000),((3000) - (100))/9)))
n_int_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)
p_X<-p_X_list[id1]
inv_r<-inv_r_list[id2]
n_int<-n_int_list[id3]
if(inv_r<10){
    stop("Finshed!!!!")
}
print(paste0("p_X=",p_X,";inv_r=",inv_r,":n_int=",n_int))
#source("/users/rzhao1/fusion/code/fusionGMM.R")
library(magic)
library(mvtnorm)
library(glmnet)
#n_int<-100
#inv_r<-3
n_ext<-n_int*inv_r
#p_X<-10
p_A<-150
XAcolnames<-c(paste0("X",1:p_X),paste0("A",1:p_A))
intercept<- -2.3
coef_X<-rep(0,p_X)
coef_A<-rep(0,p_A)
X_nonzero_index<-1:10
#A_nonzero_index<-c(11,16,17,21,26,31,34,36,50,64,67,69,88,100,200) #15 nonzero
A_nonzero_index<-c(11,12,13,17,21,26,31,34,36,50,64,67,69,88,90)-10 #15 nonzero
#coef_X[1:p_X]<-0.2
#coef_A[1:p_A]<-0.061
coef_X[X_nonzero_index]<- rnorm(length(X_nonzero_index),0.3,0.001)
coef_A[A_nonzero_index]<- rnorm(length(A_nonzero_index),0.25,0.001)
#coef_A[A_nonzero_index[1:5]]<- 0.1
#coef_A[A_nonzero_index2]<- -0.1
coefXA<-c(coef_X,coef_A)
XAc<-list()
if(F){
    XAc[[1]]<-c(8,9)
    XAc[[2]]<-c(1)
    XAc[[3]]<-c(2)
    XAc[[4]]<-c(3)
    XAc[[5]]<-c(4)
    XAc[[6]]<-c(10)
    XAc[[8]]<-c(14)
    XAc[[10]]<-c(11,13)
}else{
    XAc[[10]]<-c(8,9)
    XAc[[8]]<-c(1)
    XAc[[6]]<-c(2)
    XAc[[5]]<-c(3)
    XAc[[4]]<-c(4)
    XAc[[3]]<-c(10)
    XAc[[2]]<-c(14)
    XAc[[1]]<-c(11,13)
}
var_noise<-3
autocorM<-function(size,rho){return(rho^abs(row(diag(size))-col(diag(size))))}
generateM<-function(n,p_X,p_A,coefs,intercept,
                    rho=0.5,block_size=10,adj_r=0.3){
    #block_size=p_X+p_A
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
                #M[,jj]<-M[,jj]+adj_r*M[,i]
                M[,jj]<-M[,jj]+adj_r*M[,i]
            }
        }
    }
    #coefs<-coefs1
    #coefs1[X_nonzero_index]<-c(rep(0.1,5),rep(0.2,5))
    #M[,11]<-M%*%coefs1
    M<-scale(M)
    y<-M%*%coefs+intercept+rnorm(n,0,var_noise)

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
if(T){
    ntest = 10^4*2
    My0<-generateM(ntest,p_X,p_A,coefXA,intercept)
    M0<-My0$M
    y0<-My0$y
    #saveRDS(My0,paste0("/users/rzhao1/fusion/test/test2_",p_X,".rds"))
}else{
    My0<-readRDS(paste0("/users/rzhao1/fusion/test/test2_",p_X,".rds"))
    M0<-My0$M
    y0<-My0$y
}

R2<-function(y0,predy){
    1-mean((predy-scale(y0,scale=F))^2)/mean((scale(y0,scale = F))^2)
}

r0<-R2(y0,M0%*%coefXA)
#r0<-1-/var(y0)
r0X<-R2(y0,M0[,1:p_X]%*%coefXA[1:p_X])
r0A<-R2(y0,M0[,-c(1:p_X)]%*%coefXA[-c(1:p_X)])

My_int<-generateM(n_int,p_X,p_A,coefXA,intercept)
M_int<-My_int$M
y_int<-c(My_int$y)
y_int<-scale(y_int,scale = F)

X_int<-M_int[,1:p_X]
A_int<-M_int[,-c(1:p_X)]

fuse_lasso<-cv.glmnet(x = M_int,y = y_int)
rlasso<-R2(y0,M0%*%coef(fuse_lasso,s="lambda.min")[-1])
My_ext<-generateM(n_ext,p_X,p_A,coefXA,intercept)
M_ext<-My_ext$M
y_ext<-My_ext$y
colnames(M_ext)<-XAcolnames
sum_res<-summary_stat(M_ext,y_ext,p_X)
multilm<-sum_res$multilm
sum_mul<-sum_res$sum_mul
sum_mul[[1]]$Covariance<-sum_mul[[1]]$Covariance/sum_mul[[1]]$Sample_size
sum_uni<-sum_res$sum_uni
fuse_mul_1lam<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = T,initial_with_GMM = T)
rmul1<-R2(y0,M0%*%fuse_mul_1lam$beta)
fuse_mul_1lam_2init<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = T,initial_with_GMM = T)
rmul2<-R2(y0,M0%*%fuse_mul_1lam_2init$beta)
fuse_mul_1lam_3C<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = T,desparseC = T)
rmul3<-R2(y0,M0%*%fuse_mul_1lam_3C$beta)
fuse_mul_1lam_2init_4C<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = T,initial_with_GMM = T,desparseC = T)
rmul4<-R2(y0,M0%*%fuse_mul_1lam_2init_4C$beta)
#fuse_mul<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = T)
#rmul<-R2(y0,M0%*%fuse_mul$beta)
#c(r0,r0X,r0A,rlasso,rmul1,rmul)
print(c(r0,r0X,r0A,rlasso,rmul1,rmul2,rmul3,rmul4))
