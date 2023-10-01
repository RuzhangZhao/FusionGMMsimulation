args <- commandArgs(trailingOnly = TRUE)
id1 = as.numeric(args[1])
pacman::p_load(expm,magic,glmnet,cluster,MASS,locfit,corpcor,caret,mvtnorm)
library(intgmm)
#source("/users/rzhao1/intgmm/code/intgmm.R")
#path<-"/Users/zhaoruzhang/Desktop/UKBB_application/FusionGMMdata/nonpos/"
path<-"/users/rzhao1/intgmm/result/nonpos/"
library(magic)
library(mvtnorm)
p_X<-40
inv_r<-30
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
#ntest_list<-100*2^(1:8)*30
ntest_list<-500*2^(0:5)*30
#ntest_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)[2:10]*(1+inv_r)
ntest<-ntest_list[id1]
filename<-paste0("/users/rzhao1/intgmm/data/fdrdata/testonlyX_",p_X,"n_",ntest,".rds")
if(!file.exists(filename)){
    My0<-generateM(ntest,p_X,p_A,coefXA,intercept)
    M0<-My0$M
    saveRDS(M0,paste0("/users/rzhao1/intgmm/data/fdrdata/testonlyX_",p_X,"n_",ntest,".rds"))
}else{
    M0<-readRDS(paste0("/users/rzhao1/intgmm/data/fdrdata/testonlyX_",p_X,"n_",ntest,".rds"))
}
library(pROC)
RC<-function(pred_y){
    suppressMessages(cur_auc<-c(auc(y0,c(expit(pred_y)),direction = "<")))
}
internal_index<-c(1:floor(nrow(M0)/(1+inv_r)))
M_int<-M0[internal_index,]
M_ext<-M0[-internal_index,]
correct_pos<-c(which(coefXA!=0))

if(file.exists(paste0(path,'Lgwrap_pos_ir_',inv_r,'_pX_',p_X,'_n_',ntest,'.rds'))){
    nonpostotal<-readRDS(paste0(path,'Lgwrap_pos_ir_',inv_r,'_pX_',p_X,'_n_',ntest,'.rds'))
    APnonpostotal<-readRDS(paste0(path,'LgAPwrap_pos_ir_',inv_r,'_pX_',p_X,'_n_',ntest,'.rds'))
}else{
    nonpostotal<-list()
    APnonpostotal<-list()
}

#for(i in (length(nonpostotal)+1):100){
for(i in 1:100){
    set.seed(i)
    seed.use = sample(1:20232023,1)
    set.seed(seed.use)
    message(paste0("EPOCH",i))
    #y0<-M0%*%coefXA+intercept+rnorm(nrow(M0),0,3)
    y0<-rbinom(nrow(M0),size = 1,prob=expit(M0%*%coefXA+intercept))
    y_int<-y0[internal_index]
    y_ext<-y0[-internal_index]
    sum_res<-summary_stat(M_ext,y_ext,p_X)
    multilm<-sum_res$multilm
    sum_mul<-sum_res$sum_mul
    sum_uni<-sum_res$sum_uni

    X_int<-M_int[,1:p_X]
    A_int<-M_int[,-c(1:p_X)]
    library(intgmm)
    fuse_uni_ada<-cv.intgmm(y_int,X_int,A_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",family = "binomial",inference = T,use_sparseC = T)
    fuse_uni_ada_1lam<-cv.intgmm(y_int,X_int,A_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",family = "binomial",tune_ratio = F,inference = T,use_sparseC = T)

    fuse_mul_ada<-cv.intgmm(y_int,X_int,A_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",family = "binomial",inference = T,use_sparseC = T)
    fuse_mul_ada_1lam<-cv.intgmm(y_int,X_int,A_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",family = "binomial",tune_ratio = F,inference = T,use_sparseC = T)

    y_int_g<-scale(y_int,scale = F)/var(y_int)
    APfuse_uni_ada<-cv.intgmm(y_int_g,X_int,A_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",inference = T,use_sparseC = T)
    APfuse_uni_ada_1lam<-cv.intgmm(y_int_g,X_int,A_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",tune_ratio = F,inference = T,use_sparseC = T)

    APfuse_mul_ada<-cv.intgmm(y_int_g,X_int,A_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",inference = T,use_sparseC = T)
    APfuse_mul_ada_1lam<-cv.intgmm(y_int_g,X_int,A_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",tune_ratio = F,inference = T,use_sparseC = T)


    #fuse_ridge<-glm(y~.,data.frame(y=y_int,M_int),family = "binomial")
    fuse_ridge<-cv.glmnet(x = M_int,y = y_int,family = "binomial",alpha=0)
    #www<-1/abs(c(fuse_ridge$coefficients[-1]))^(1/2)
    www<-1/abs(c(coef(fuse_ridge,s="lambda.min")[-1]))^(1/2)
    www[is.infinite(www)]<-max(www[!is.infinite(www)])*100
    fuse_ada<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www,family = "binomial")

    ada_res_func<-function(fuse_ada){
        beta_ada<-as.vector(coef(fuse_ada,s="lambda.min"))
        beta_ada<-c(beta_ada[-1],beta_ada[1])

        diag_term<-c(expit(cbind(M_int,1)%*%beta_ada)*(1-expit(cbind(M_int,1)%*%beta_ada)))
        Sigsum_scaled<-t(cbind(M_int,1))%*%(diag_term*cbind(M_int,1))#/nrow(M_int)
        index_nonzero<-which(beta_ada!=0)
        Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
        inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
        final_v<-diag(inv_Sigsum_scaled_nonzero)
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

    fdrgmm<-function(fuse_uni_ada){
        pval<-fuse_uni_ada$pval
        cp0<-fuse_uni_ada$nonzero_pos
        pval<-pval[cp0%in%c(1:(p_X+p_A))]
        cp0<-cp0[cp0%in%c(1:(p_X+p_A))]
        pval1<-p.adjust(pval,method = "BH")
        sum(!cp0[pval1<0.05]%in%correct_pos)/max(1e-8,length(cp0[pval1<0.05]))
    }
    fdrgmm0<-function(fuse_uni_ada){
        cp0<-fuse_uni_ada[['corrected_pos']]
        sum(!cp0%in%correct_pos)/max(1e-8,length(cp0))
    }
    power<-function(fuse_uni_ada){
        pval<-fuse_uni_ada$pval
        cp0<-fuse_uni_ada$nonzero_pos
        pval<-pval[cp0%in%c(1:(p_X+p_A))]
        cp0<-cp0[cp0%in%c(1:(p_X+p_A))]
        pval1<-p.adjust(pval,method = "BH")
        sum(cp0[pval1<0.05]%in%correct_pos)/length(correct_pos)
        #sum(fuse_uni_ada[['corrected_pos']]%in%correct_pos)/length(correct_pos)
    }
    power0<-function(fuse_uni_ada){
        cp0<-fuse_uni_ada[['nonzero_pos']][fuse_uni_ada[['pval']]<0.05]
        sum(cp0%in%correct_pos)/length(correct_pos)
    }
    cover1<-function(fuse_uni_ada){
        coefXA1<-c(coefXA,intercept)
        ans<-rep(F,length(correct_pos))
        for(i in 1:length(correct_pos)){
            pos = correct_pos[i]
            if( fuse_uni_ada$beta[pos]!=0){
                nonzerovar<-fuse_uni_ada$nonzero_var
                names(nonzerovar)<-fuse_uni_ada$nonzero_pos
                up<-fuse_uni_ada$beta[pos]+1.96*sqrt(nonzerovar[paste0(pos)])
                low<-fuse_uni_ada$beta[pos]-1.96*sqrt(nonzerovar[paste0(pos)])
                if( coefXA1[pos]<up &coefXA1[pos]>low  ){
                    ans[i]<-T
                }else{
                    print("Notcover!")
                }
            }
        }
        mean(ans)
    }
    cover<-function(fuse_uni_ada){
        #coefXA1<-c(coefXA,intercept)
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

    cover2<-function(fuse_uni_ada){
        #coefXA1<-c(coefXA,intercept)
        ans<-rep(F,length(correct_pos))
        for(i in 1:length(correct_pos)){
            pos = correct_pos[i]
            if( fuse_uni_ada$beta[pos]!=0){
                nonzerovar<-fuse_uni_ada$nonzero_var
                names(nonzerovar)<-fuse_uni_ada$nonzero_pos
                up<-fuse_uni_ada$beta[pos]+2*sqrt(nonzerovar[paste0(pos)])
                low<-fuse_uni_ada$beta[pos]-2*sqrt(nonzerovar[paste0(pos)])
                if( coefXA[pos]<up &coefXA[pos]>low  ){
                    ans[i]<-T
                }else{
                    print("Notcover!")
                }
            }
        }
        mean(ans)
    }



    pvalfunc<-function(ada_res,fuse_uni_ada,fuse_uni_ada_1lam,fuse_mul_ada,fuse_mul_ada_1lam){
        nonpos<-c(fdrgmm(ada_res),fdrgmm(fuse_uni_ada),fdrgmm(fuse_uni_ada_1lam),fdrgmm(fuse_mul_ada),fdrgmm(fuse_mul_ada_1lam))
        #nonpos0<-c(fdrgmm0(ada_res),fdrgmm0(fuse_uni_ada),fdrgmm0(fuse_uni_ada_1lam),fdrgmm0(fuse_mul_ada),fdrgmm0(fuse_mul_ada_1lam))
        powers<-c(power(ada_res),power(fuse_uni_ada),power(fuse_uni_ada_1lam),power(fuse_mul_ada),power(fuse_mul_ada_1lam))
        #powers0<-c(power0(ada_res),power0(fuse_uni_ada),power0(fuse_uni_ada_1lam),power0(fuse_mul_ada),power0(fuse_mul_ada_1lam))
        #cover2age<-c(cover2(ada_res),cover2(fuse_uni_ada),cover2(fuse_uni_ada_1lam),cover2(fuse_mul_ada),cover2(fuse_mul_ada_1lam))
        #cover1age<-c(cover1(ada_res),cover1(fuse_uni_ada),cover1(fuse_uni_ada_1lam),cover1(fuse_mul_ada),cover1(fuse_mul_ada_1lam))
        coverage<-c(cover(ada_res),cover(fuse_uni_ada),cover(fuse_uni_ada_1lam),cover(fuse_mul_ada),cover(fuse_mul_ada_1lam))
        allname<-c("ada","uni_ada","uni_ada_1lam","mul_ada","mul_ada_1lam")

        pval_list<-list()
        pval_list[[1]]<-ada_res$pval
        pval_list[[2]]<-fuse_uni_ada$pval
        pval_list[[3]]<-fuse_uni_ada_1lam$pval
        pval_list[[4]]<-fuse_mul_ada$pval
        pval_list[[5]]<-fuse_mul_ada_1lam$pval
        nonzero_pos_list<-list()
        nonzero_pos_list[[1]]<-ada_res$nonzero_pos
        nonzero_pos_list[[2]]<-fuse_uni_ada$nonzero_pos
        nonzero_pos_list[[3]]<-fuse_uni_ada_1lam$nonzero_pos
        nonzero_pos_list[[4]]<-fuse_mul_ada$nonzero_pos
        nonzero_pos_list[[5]]<-fuse_mul_ada_1lam$nonzero_pos
        names(nonpos)<-allname
        #names(nonpos0)<-allname
        names(powers)<-allname
        #names(powers0)<-allname
        names(coverage)<-allname
        #names(cover1age)<-allname
        #names(cover2age)<-allname
        print("FDR")
        print(c(nonpos))
        print("Power")
        print(c(powers))
        print("Cover")
        print(c(coverage))
        #list(nonpos,nonpos0,powers,powers0,coverage,pval_list,nonzero_pos_list)
        list(nonpos,powers,coverage,pval_list,nonzero_pos_list)
    }


    nonpostotal[[i]]<-pvalfunc(ada_res,fuse_uni_ada,fuse_uni_ada_1lam,fuse_mul_ada,fuse_mul_ada_1lam)
    APnonpostotal[[i]]<-pvalfunc(ada_res,APfuse_uni_ada,APfuse_uni_ada_1lam,APfuse_mul_ada,APfuse_mul_ada_1lam)

    if (i %% 5 == 0){
        saveRDS(nonpostotal,paste0(path,'Lgwrap_pos_ir_',inv_r,'_pX_',p_X,'_n_',ntest,'.rds'))
        saveRDS(APnonpostotal,paste0(path,'LgAPwrap_pos_ir_',inv_r,'_pX_',p_X,'_n_',ntest,'.rds'))
    }
}



