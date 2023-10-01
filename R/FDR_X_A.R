source("/users/rzhao1/intgmm/code/intgmm.R")
#path<-"/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/nonpos/"
path<-"/users/rzhao1/intgmm/result/nonpos/"
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
    ntest = 40000
    My0<-generateM(ntest,p_X,p_A,coefXA,intercept)
    M0<-My0$M
    saveRDS(M0,paste0("/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/testonlyX_",p_X,".rds"))
}else{
    #M0<-readRDS(paste0("/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/testonlyX_",p_X,".rds"))
    M0<-readRDS(paste0("/users/rzhao1/intgmm/data/fdrdata/testonlyX_",p_X,".rds"))
}

R2<-function(y0,predy){
    1-sum((predy-scale(y0,scale = F))^2)/sum((scale(y0,scale = F))^2)
}
internal_index<-c(1:floor(nrow(M0)/(1+inv_r)))
M_int<-M0[internal_index,]
M_ext<-M0[-internal_index,]
correct_pos<-which(coefXA!=0)
nonpostotal<-list()
for(i in 1:100){
    seed.use = sample(1:20232023,1)
    set.seed(seed.use)
    message(paste0("EPOCH",i))
    y0<-M0%*%coefXA+intercept+rnorm(nrow(M0),0,3)
    R2(y0,M0%*%coefXA)
    y_int<-y0[internal_index]
    y_ext<-y0[-internal_index]
    sum_res<-summary_stat(M_ext,y_ext,p_X)
    multilm<-sum_res$multilm
    sum_mul<-sum_res$sum_mul
    sum_uni<-sum_res$sum_uni

    X_int<-M_int[,1:p_X]
    A_int<-M_int[,-c(1:p_X)]

    fuse_uni_ada<-fusionGMM.addition(X_int,A_int,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",inference = T,seed.use=seed.use)
    fuse_uni_ada_1lam<-fusionGMM.addition(X_int,A_int,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",tune_ratio = F,inference = T,seed.use=seed.use)

    fuse_mul_ada<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",inference = T,seed.use=seed.use)
    fuse_mul_ada_1lam<-fusionGMM.addition(X_int,A_int,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",tune_ratio = F,inference = T,seed.use=seed.use)


    fuse_ridge<-lm(y~.,data.frame(y=y_int,M_int))
        #cv.glmnet(x = M_int,y = y_int,alpha=0)
    www<-1/abs(c(fuse_ridge$coefficients[-1]))^(1/2)
    www1<-1/abs(c(fuse_ridge$coefficients[-1]))
    www2<-1/abs(c(fuse_ridge$coefficients[-1]))^2
    #www<-1/abs(c(coef(fuse_ridge,s="lambda.min")[-1]))^(1/2)
    #www[is.infinite(www)]<-max(www[!is.infinite(www)])*100
    fuse_ada<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www)
    fuse_ada1<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www1)
    fuse_ada2<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www2)

    ada_res_func<-function(fuse_ada){
        beta_ada<-coef(fuse_ada,s="lambda.min")[-1]
        sigma_hat<-mean((M_int%*%beta_ada-y_int)^2)
        index_nonzero<-which(beta_ada!=0)
        Sigsum_scaled_nonzero<-crossprod(M_int[,index_nonzero])/nrow(M_int)
        inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
        final_v<-diag(inv_Sigsum_scaled_nonzero)*sigma_hat

        pval_final<-pchisq(nrow(M_int)*beta_ada[index_nonzero]^2/final_v,1,lower.tail = F)
        #corrected_pos<-index_nonzero[which(pval_final<0.05/length(index_nonzero))]
        pval_final1<-p.adjust(pval_final,method = "BH")
        corrected_pos<-index_nonzero[which(pval_final1<0.05)]
        corrected_pos0<-index_nonzero[which(pval_final<0.05)]
        ada_res<-list()
        ada_res$corrected_pos<-corrected_pos
        ada_res$pval<-pval_final
        ada_res$beta<-beta_ada
        ada_res$nonzero_var<-final_v
        ada_res$nonzero_pos<-index_nonzero
        ada_res
    }
    ada_res<-ada_res_func(fuse_ada)
    ada_res1<-ada_res_func(fuse_ada1)
    ada_res2<-ada_res_func(fuse_ada2)

    fdrgmm<-function(fuse_uni_ada){
        sum(!fuse_uni_ada[['corrected_pos']]%in%correct_pos)/length(fuse_uni_ada[['corrected_pos']])
    }
    fdrgmm0<-function(fuse_uni_ada){
        cp0<-fuse_uni_ada[['nonzero_pos']][fuse_uni_ada[['pval']]<0.05]
        sum(!cp0%in%correct_pos)/length(cp0)
    }
    power<-function(fuse_uni_ada){
        sum(fuse_uni_ada[['corrected_pos']]%in%correct_pos)/length(correct_pos)
    }
    power0<-function(fuse_uni_ada){
        cp0<-fuse_uni_ada[['nonzero_pos']][fuse_uni_ada[['pval']]<0.05]
        sum(cp0%in%correct_pos)/length(correct_pos)
    }
    cover<-function(fuse_uni_ada){
        ans<-rep(F,length(correct_pos))
        for(i in 1:length(correct_pos)){
            pos = correct_pos[i]
            if( fuse_uni_ada$beta[pos]!=0){
                nonzerovar<-fuse_uni_ada$nonzero_var
                names(nonzerovar)<-fuse_uni_ada$nonzero_pos
                up<-fuse_uni_ada$beta[pos]+1.96*sqrt(nonzerovar[paste0(pos)])
                low<-fuse_uni_ada$beta[pos]-1.96*sqrt(nonzerovar[paste0(pos)])
                if( coefXA[pos]<up &coefXA[pos]>low  ){
                    ans[i]<-T
                }else{
                    print("Notcover!")
                }
            }
        }
        mean(ans)
    }

    nonpos<-c(fdrgmm(ada_res),fdrgmm(ada_res1),fdrgmm(ada_res2),fdrgmm(fuse_uni_ada),fdrgmm(fuse_uni_ada_1lam),fdrgmm(fuse_mul_ada),fdrgmm(fuse_mul_ada_1lam))
    nonpos0<-c(fdrgmm0(ada_res),fdrgmm0(ada_res1),fdrgmm0(ada_res2),fdrgmm0(fuse_uni_ada),fdrgmm0(fuse_uni_ada_1lam),fdrgmm0(fuse_mul_ada),fdrgmm0(fuse_mul_ada_1lam))
    powers<-c(power(ada_res),power(ada_res1),power(ada_res2),power(fuse_uni_ada),power(fuse_uni_ada_1lam),power(fuse_mul_ada),power(fuse_mul_ada_1lam))
    powers0<-c(power0(ada_res),power0(ada_res1),power0(ada_res2),power0(fuse_uni_ada),power0(fuse_uni_ada_1lam),power0(fuse_mul_ada),power0(fuse_mul_ada_1lam))
    coverage<-c(cover(ada_res),cover(ada_res1),cover(ada_res2),cover(fuse_uni_ada),cover(fuse_uni_ada_1lam),cover(fuse_mul_ada),cover(fuse_mul_ada_1lam))
    allname<-c("ada","ada1","ada2","uni_ada","uni_ada_1lam","mul_ada","mul_ada_1lam")
    names(nonpos)<-allname
    names(nonpos0)<-allname
    names(powers)<-allname
    names(powers0)<-allname
    names(coverage)<-allname
    nonpostotal[[i]]<-list(nonpos,nonpos0,powers,powers0,coverage)
    print(c(nonpos,nonpos0))
    print(c(powers,powers0))
    print(coverage)
    if (i %% 5 == 0){
        saveRDS(nonpostotal,paste0(path,'pos_ir_',inv_r,'_pX_',p_X,'2.rds'))
    }
}



