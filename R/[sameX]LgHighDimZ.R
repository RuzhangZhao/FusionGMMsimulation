args <- commandArgs(trailingOnly = TRUE)
id3s = as.numeric(args[1])
use_sparseC=TRUE
library(htlgmm)
pacman::p_load(expm,magic,glmnet,cluster,MASS,locfit,corpcor,caret,pROC,mvtnorm)
p_X_list<-c(10,40,150)
n_int_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)
foldname="_same"
filecount=""

for(id3 in c(id3s)){
    p_X<-150
    n_int<-n_int_list[id3]
    print(paste0("p_X=",p_X,";inv_r=",10,30,":n_int=",n_int))
    #source("/users/rzhao1/fusion/code/fusionGMM.R")
    library(magic)
    library(mvtnorm)
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
        y<-rbinom(n,size = 1,prob=expit(M%*%coefs+intercept))
        list("M"=M,"y"=y)
    }


    ada_res_func<-function(fuse_ada,M_int){
        beta_ada<-as.vector(coef(fuse_ada,s="lambda.min"))
        beta_ada<-c(beta_ada[-1],beta_ada[1])

        diag_term<-c(expit(cbind(M_int,1)%*%beta_ada)*(1-expit(cbind(M_int,1)%*%beta_ada)))
        Sigsum_scaled<-t(cbind(M_int,1))%*%(diag_term*cbind(M_int,1))#/nrow(M_int)
        index_nonzero<-which(beta_ada!=0)
        Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
        final_v<-solve(Sigsum_scaled_nonzero)

        ada_res<-list()
        ada_res$beta<-beta_ada[-length(beta_ada)]
        ada_res$var<-final_v[-length(index_nonzero),-length(index_nonzero)]
        ada_res$nonzero_pos<-index_nonzero[-length(index_nonzero)]
        ada_res
    }

    summary_stat<-function(M,y,p_X){
        study_info_scaled<-list()
        library(glmnet)
        fuse_ridge<-cv.glmnet(x = M[,1:p_X],y = y,family = "binomial",alpha=0)
        www<-1/abs(c(coef(fuse_ridge,s="lambda.min")[-1]))^(1/2)
        www[is.infinite(www)]<-max(www[!is.infinite(www)])*100
        fuse_ada<-cv.glmnet(x = M[,1:p_X],y = y,alpha=1,penalty.factor = www,family = "binomial")

        ada_res<-ada_res_func(fuse_ada,M[,1:p_X])
        ada_var<-ada_res$var
        nonzero_pos = ada_res$nonzero_pos
        study.m = list(Coeff=ada_res$beta[ada_res$nonzero_pos],
                       Covariance=ada_var,
                       Sample_size=nrow(M))
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
        return(list("multilm"=fuse_ada,"sum_mul"=study_info_scaled,"sum_uni"=study_info_scaled_sep,"nonzero_pos"=nonzero_pos))
    }

    ### Infinity Test Data for Validation
    if(F){
        ntest = 10^6
        #My0<-generateM(ntest,p_X,p_A,coefXA,intercept)
        #saveRDS(My0,paste0("/users/rzhao1/fusion/test/test2_",p_X,".rds"))
        #saveRDS(My0,paste0("../FusionGMMdata/test2_",p_X,".rds"))
        My0<-readRDS(paste0("../FusionGMMdata/test2_",p_X,".rds"))
        M0<-My0$M
        y0=rbinom(nrow(M0),size = 1,prob=expit(M0%*%coefXA+intercept))
        saveRDS(y0,paste0("/users/rzhao1/intgmm/data/test2binary_y_",p_X,".rds"))
        saveRDS(y0,paste0("../FusionGMMdata/test2binary_y_",p_X,".rds"))
    }else{
        #My0<-readRDS(paste0("/users/rzhao1/fusion/test/test2_",p_X,".rds"))
        My0<-readRDS(paste0("/users/rzhao1/intgmm/data/test2_",p_X,".rds"))
        M0<-My0$M
        y0<-readRDS(paste0("/users/rzhao1/intgmm/data/test2binary_y_",p_X,".rds"))
    }

    R2<-function(y0,predy){
        1-sum((predy-scale(y0,scale = F))^2)/sum((scale(y0,scale = F))^2)
    }
    if(!file.exists(paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/Lg_n_',n_int,'_ir_',10,'_pX_',p_X,filecount,'.rds'))){
        rr_gap<-list()
        rr_gaplinear<-list()
        rr_gap2<-list()
        rr_gap2linear<-list()
    }else{
        rr_gap<-readRDS(paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/Lg_n_',n_int,'_ir_',10,'_pX_',p_X,filecount,'.rds'))
        #rr_gaplinear<-readRDS(paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/LgAP_n_',n_int,'_ir_',10,'_pX_',p_X,filecount,'.rds'))
        rr_gap2<-readRDS(paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/Lg_n_',n_int,'_ir_',30,'_pX_',p_X,filecount,'.rds'))
        #rr_gap2linear<-readRDS(paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/LgAP_n_',n_int,'_ir_',30,'_pX_',p_X,filecount,'.rds'))
    }
    if(length(rr_gap)==100){break}
    for(i in (length(rr_gap)+1):100){
        set.seed(i)
        seed.use<-sample(1:20232023,1)
        set.seed(seed.use)
        if(filecount=='_1'){
            seed.use<-sample(1:2147483647,1)
            set.seed(seed.use)
        }
        if(filecount=='_2'){
            seed.use<-sample(1:20232023,1)
            set.seed(seed.use)
            seed.use<-sample(1:2147483647,1)
            set.seed(seed.use)
        }
        if(filecount=='_3'){
            seed.use<-sample(1:20232023,1)
            set.seed(seed.use)
            seed.use<-sample(1:20232023,1)
            set.seed(seed.use)
            seed.use<-sample(1:2147483647,1)
            set.seed(seed.use)
        }
        message(paste0("EPOCH",i,"seed:",seed.use))
        # external data
        My_ext<-generateM(n_int*10,p_X,p_A,coefXA,intercept)
        M_ext<-My_ext$M
        y_ext<-My_ext$y
        colnames(M_ext)<-XAcolnames
        sum_res<-summary_stat(M_ext,y_ext,p_X)
        multilm<-sum_res$multilm
        sum_mul<-sum_res$sum_mul
        sum_uni<-sum_res$sum_uni
        nonzero_pos=sum_res$nonzero_pos

        My_ext<-generateM(n_int*30,p_X,p_A,coefXA,intercept)
        M_ext<-My_ext$M
        y_ext<-My_ext$y
        colnames(M_ext)<-XAcolnames
        sum_res<-summary_stat(M_ext,y_ext,p_X)
        multilm2<-sum_res$multilm
        sum_mul2<-sum_res$sum_mul
        sum_uni2<-sum_res$sum_uni
        nonzero_pos2=sum_res$nonzero_pos
        # internal data
        My_int<-generateM(n_int,p_X,p_A,coefXA,intercept)
        M_int<-My_int$M
        y_int<-c(My_int$y)
        X_int<-M_int[,1:p_X]
        A_int<-M_int[,-c(1:p_X)]
        y_int_g<-scale(y_int,scale = F)/var(y_int)

        fuse_lasso<-cv.glmnet(x = M_int,y = y_int,family = "binomial")
        fuse_ridge<-cv.glmnet(x = M_int,y = y_int,alpha=0,family = "binomial")
        www<-1/abs(c(coef(fuse_ridge,s="lambda.min")[-1]))^(1/2)
        www[is.infinite(www)]<-max(www[!is.infinite(www)])*100
        fuse_ada<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www,family = "binomial")
        cat("5%")
        library(pROC)
        RC<-function(pred_y){
            suppressMessages(cur_auc<-c(auc(y0,c(expit(pred_y)),direction = "<")))
        }

        r0<-RC(M0%*%coefXA)
        r0X<-RC(M0[,1:p_X]%*%coef_X)
        r0A<-RC(M0[,-c(1:p_X)]%*%coef_A)
        rlasso<-RC(M0%*%coef(fuse_lasso,s="lambda.min")[-1])
        rlasso1<-RC(M0%*%coef(fuse_lasso,s="lambda.1se")[-1])
        rada<-RC(M0%*%coef(fuse_ada,s="lambda.min")[-1])
        rada1<-RC(M0%*%coef(fuse_ada,s="lambda.1se")[-1])
        runhtlgmm<-function(y_int,X_int,A_int,sum_uni,sum_mul,family="binomial",nonzero_pos){
            #fuse_uni<-cv.htlgmm(y_int,X_int,A_int,sum_uni,penalty_type = "lasso",summary_type = "uni",family=family,tune_ratio = T,use_sparseC = use_sparseC)
            fuse_uni_1lam<-cv.htlgmm(y_int,X_int,A_int,sum_uni,penalty_type = "lasso",summary_type = "uni",family=family,tune_ratio = F,use_sparseC = use_sparseC)
            fuse_uni=fuse_uni_1lam
            cat("20%")
            #fuse_uni_ada<-cv.htlgmm(y_int,X_int,A_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",family=family,tune_ratio = T,use_sparseC = use_sparseC)
            fuse_uni_ada_1lam<-cv.htlgmm(y_int,X_int,A_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",family=family,tune_ratio = F,use_sparseC = use_sparseC)
            fuse_uni_ada=fuse_uni_ada_1lam
            cat("40%")
            fuse_mul<-cv.htlgmm(y_int,X_int[,nonzero_pos],A_int,sum_mul,penalty_type = "lasso",summary_type = "multi",family=family,tune_ratio = T,use_sparseC = use_sparseC)
            fuse_mul_1lam<-cv.htlgmm(y_int,X_int[,nonzero_pos],A_int,sum_mul,penalty_type = "lasso",summary_type = "multi",family=family,tune_ratio = F,use_sparseC = use_sparseC)
            cat("60%")
            fuse_mul_ada<-cv.htlgmm(y_int,X_int[,nonzero_pos],A_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",family=family,tune_ratio = T,use_sparseC = use_sparseC)
            fuse_mul_ada_1lam<-cv.htlgmm(y_int,X_int[,nonzero_pos],A_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",family=family,tune_ratio = F,use_sparseC = use_sparseC)
            cat("80%")
            rX<-RC(M0[,1:p_X][,nonzero_pos]%*%sum_mul[[1]]$Coeff)
            runi<-RC(M0%*%fuse_uni$beta[1:(p_X+p_A)])
            runi1<-RC(M0%*%fuse_uni_1lam$beta[1:(p_X+p_A)])
            runiada<-RC(M0%*%fuse_uni_ada$beta[1:(p_X+p_A)])
            runiada1<-RC(M0%*%fuse_uni_ada_1lam$beta[1:(p_X+p_A)])

            rmul<-RC(M0[,c(nonzero_pos,c((p_X+1):(p_X+p_A)))]%*%fuse_mul$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))])
            rmul1<-RC(M0[,c(nonzero_pos,c((p_X+1):(p_X+p_A)))]%*%fuse_mul_1lam$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))])
            rmulada<-RC(M0[,c(nonzero_pos,c((p_X+1):(p_X+p_A)))]%*%fuse_mul_ada$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))])
            rmulada1<-RC(M0[,c(nonzero_pos,c((p_X+1):(p_X+p_A)))]%*%fuse_mul_ada_1lam$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))])
            betagap<-c(sum((coef(fuse_lasso,s="lambda.min")[-1] - coefXA)^2),
                       sum((coef(fuse_lasso,s="lambda.1se")[-1] - coefXA)^2),
                       sum((coef(fuse_ada,s="lambda.min")[-1] - coefXA)^2),
                       sum((coef(fuse_ada,s="lambda.1se")[-1] - coefXA)^2),
                       sum((fuse_mul$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))] - coefXA[c(nonzero_pos,c((p_X+1):(p_X+p_A)))])^2),
                       sum((fuse_mul_1lam$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))] - coefXA[c(nonzero_pos,c((p_X+1):(p_X+p_A)))])^2),
                       sum((fuse_mul_ada$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))] - coefXA[c(nonzero_pos,c((p_X+1):(p_X+p_A)))])^2),
                       sum((fuse_mul_ada_1lam$beta[1:(length(c(nonzero_pos,c((p_X+1):(p_X+p_A)))))] - coefXA[c(nonzero_pos,c((p_X+1):(p_X+p_A)))])^2),
                       sum((fuse_uni$beta[1:(p_X+p_A)] - coefXA)^2),
                       sum((fuse_uni_1lam$beta[1:(p_X+p_A)] - coefXA)^2),
                       sum((fuse_uni_ada$beta[1:(p_X+p_A)] - coefXA)^2),
                       sum((fuse_uni_ada_1lam$beta[1:(p_X+p_A)] - coefXA)^2)
            )
            names(betagap)<-c("rlasso","rlasso1","rada","rada1","rmul","rmul1","rmulada","rmulada1","runi","runi1","runiada","runiada1")
            cat("curres\n")
            rr<-c(r0,r0X,r0A,rX,rlasso,rlasso1,rada,rada1,rmul,rmul1,rmulada,rmulada1,runi,runi1,runiada,runiada1)
            names(rr)<-c("r0","r0X","r0A","rX","rlasso","rlasso1","rada","rada1","rmul","rmul1","rmulada","rmulada1","runi","runi1","runiada","runiada1")
            print(round(rr,3))
            print(round(betagap,3))
            ratios<-c(fuse_uni$ratio_min,fuse_uni_ada$ratio_min,fuse_mul$ratio_min,
                      fuse_mul_ada$ratio_min)
            print(round(ratios,3))
            list("rr"=rr,"beta"=betagap,"ratios"=ratios)
        }
        print("inv_r=10+Logistic Quadratic")
        rr_gap[[i]]<-runhtlgmm(y_int,X_int,A_int,sum_uni,sum_mul,family="binomial",nonzero_pos)
        #print("inv_r=10+Logistic Linear")
        #rr_gaplinear[[i]]<-runhtlgmm(y_int_g,X_int,A_int,sum_uni,sum_mul,family="gaussian",nonzero_pos)
        print("inv_r=30+Logistic Quadratic")
        rr_gap2[[i]]<-runhtlgmm(y_int,X_int,A_int,sum_uni2,sum_mul2,family="binomial",nonzero_pos2)
        #print("inv_r=30+Logistic Linear")
        #rr_gap2linear[[i]]<-runhtlgmm(y_int_g,X_int,A_int,sum_uni2,sum_mul2,family="gaussian",nonzero_pos)

        cat("100%\n")
        if (i %% 5 == 0){
            saveRDS(rr_gap,paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/Lg_n_',n_int,'_ir_',10,'_pX_',p_X,filecount,'.rds'))
            #saveRDS(rr_gaplinear,paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/LgAP_n_',n_int,'_ir_',10,'_pX_',p_X,filecount,'.rds'))
            saveRDS(rr_gap2,paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/Lg_n_',n_int,'_ir_',30,'_pX_',p_X,filecount,'.rds'))
            #saveRDS(rr_gap2linear,paste0('/users/rzhao1/intgmm/result',foldname,'/Lg/LgAP_n_',n_int,'_ir_',30,'_pX_',p_X,filecount,'.rds'))
        }
    }
}
