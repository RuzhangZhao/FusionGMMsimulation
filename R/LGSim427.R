library(locfit)
library(magic)
library(mvtnorm)
for(id2 in c(1)){
    for(id1 in c(2)){
        for(id3 in c(3:10)){
            p_X_list<-c(10,40)
            inv_r_list<-c(10,30)
            n_int_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)
            p_X<-p_X_list[id1]
            inv_r<-inv_r_list[id2]
            n_int<-n_int_list[id3]
            print(paste0("p_X=",p_X,";inv_r=",inv_r,":n_int=",n_int))
            #source("/users/rzhao1/fusion/code/fusionGMM.R")
            library(magic)
            library(mvtnorm)
            #n_int<-100
            #inv_r<-3
            n_ext<-n_int*inv_r
            #p_X<-10
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
            if(F){
                ntest = 10^6
                #My0<-generateM(ntest,p_X,p_A,coefXA,intercept)
                #saveRDS(My0,paste0("/users/rzhao1/fusion/test/test2_",p_X,".rds"))
                #saveRDS(My0,paste0("../FusionGMMdata/test2_",p_X,".rds"))
                My0<-readRDS(paste0("../FusionGMMdata/test2_",p_X,".rds"))
                M0<-My0$M
                y0=rbinom(nrow(M0),size = 1,prob=expit(M0%*%coefXA+intercept))
                saveRDS(y0,paste0("../FusionGMMdata/test2binary_y_",p_X,".rds"))
            }else{
                #My0<-readRDS(paste0("/users/rzhao1/fusion/test/test2_",p_X,".rds"))
                My0<-readRDS(paste0("../FusionGMMdata/test2_",p_X,".rds"))
                M0<-My0$M
                y0<-readRDS(paste0("../FusionGMMdata/test2binary_y_",p_X,".rds"))
            }

            R2<-function(y0,predy){
                1-sum((predy-scale(y0,scale = F))^2)/sum((scale(y0,scale = F))^2)
            }
            rr_gap<-list()
            for(i in 1:10){
                set.seed(sample(1:2023,1))
                message(paste0("EPOCH",i))
                # external data
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
                X_int<-M_int[,1:p_X]
                A_int<-M_int[,-c(1:p_X)]

                fuse_uni<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_uni,penalty_type = "lasso",summary_type = "uni",approx_cross_validation =F)
                fuse_uni_1lam<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_uni,penalty_type = "lasso",summary_type = "uni",approx_cross_validation =F,tune_ratio = F)
                cat("20%")
                fuse_uni_ada<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",approx_cross_validation =F)
                fuse_uni_ada_1lam<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_uni,penalty_type = "adaptivelasso",summary_type = "uni",approx_cross_validation =F,tune_ratio = F)
                cat("40%")
                fuse_mul<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F)
                fuse_mul_1lam<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_mul,penalty_type = "lasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = F)
                cat("60%")
                fuse_mul_ada<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",approx_cross_validation =F)
                fuse_mul_ada_1lam<-fusionGMM.binary.addition(X_int,A_int,1,y_int,sum_mul,penalty_type = "adaptivelasso",summary_type = "multi",approx_cross_validation =F,tune_ratio = F)
                cat("80%")

                fuse_lasso<-cv.glmnet(x = M_int,y = y_int,family = "binomial")
                fuse_ridge<-cv.glmnet(x = M_int,y = y_int,alpha=0,family = "binomial")
                www<-1/abs(c(coef(fuse_ridge,s="lambda.min")[-1]))^(1/2)
                www[is.infinite(www)]<-max(www[!is.infinite(www)])*100
                fuse_ada<-cv.glmnet(x = M_int,y = y_int,penalty.factor = www,family = "binomial")
                cat("95%\n")
                library(pROC)
                RC<-function(pred_y){
                    suppressMessages(cur_auc<-c(auc(y0,expit(pred_y),direction = "<")))
                }

                r0<-RC(M0%*%coefXA)
                r0X<-RC(M0[,1:p_X]%*%coef_X)
                r0A<-RC(M0[,-c(1:p_X)]%*%coef_A)
                rX<-RC(M0[,1:p_X]%*%sum_mul[[1]]$Coeff)
                rlasso<-RC(M0%*%coef(fuse_lasso,s="lambda.min")[-1])
                rlasso1<-RC(M0%*%coef(fuse_lasso,s="lambda.1se")[-1])
                rada<-RC(M0%*%coef(fuse_ada,s="lambda.min")[-1])
                rada1<-RC(M0%*%coef(fuse_ada,s="lambda.1se")[-1])

                runi<-RC(M0%*%fuse_uni$beta)
                runi1<-RC(M0%*%fuse_uni_1lam$beta)
                runiada<-RC(M0%*%fuse_uni_ada$beta)
                runiada1<-RC(M0%*%fuse_uni_ada_1lam$beta)

                rmul<-RC(M0%*%fuse_mul$beta)
                rmul1<-RC(M0%*%fuse_mul_1lam$beta)
                rmulada<-RC(M0%*%fuse_mul_ada$beta)
                rmulada1<-RC(M0%*%fuse_mul_ada_1lam$beta)

                #rmulap<-R2(y0,M0%*%fuse_mul_app$beta)
                #rmul1ap<-R2(y0,M0%*%fuse_mul_1lam_app$beta)

                betagap<-c(sum((coef(fuse_lasso,s="lambda.min")[-1] - coefXA)^2),
                           sum((coef(fuse_lasso,s="lambda.1se")[-1] - coefXA)^2),
                           sum((coef(fuse_ada,s="lambda.min")[-1] - coefXA)^2),
                           sum((coef(fuse_ada,s="lambda.1se")[-1] - coefXA)^2),
                           sum((fuse_mul$beta - coefXA)^2),
                           sum((fuse_mul_1lam$beta - coefXA)^2),
                           sum((fuse_mul_ada$beta - coefXA)^2),
                           sum((fuse_mul_ada_1lam$beta - coefXA)^2),
                           sum((fuse_uni$beta - coefXA)^2),
                           sum((fuse_uni_1lam$beta - coefXA)^2),
                           sum((fuse_uni_ada$beta - coefXA)^2),
                           sum((fuse_uni_ada_1lam$beta - coefXA)^2)
                )
                names(betagap)<-c("rlasso","rlasso1","rada","rada1","rmul","rmul1","rmulada","rmulada1","runi","runi1","runiada","runiada1")

                rr<-c(r0,r0X,r0A,rX,rlasso,rlasso1,rada,rada1,rmul,rmul1,rmulada,rmulada1,runi,runi1,runiada,runiada1)
                names(rr)<-c("r0","r0X","r0A","rX","rlasso","rlasso1","rada","rada1","rmul","rmul1","rmulada","rmulada1","runi","runi1","runiada","runiada1")
                print(round(rr,3))
                print(round(betagap,3))
                ratios<-c(fuse_uni$ratio_min,fuse_uni_ada$ratio_min,fuse_mul$ratio_min,
                          fuse_mul_ada$ratio_min)
                print(round(ratios,3))
                rr_gap[[i]]<-list("rr"=rr,"beta"=betagap,"ratios"=ratios)

                if (i %% 5 == 0){
                    saveRDS(rr_gap,paste0('../FusionGMMdata/Lg/Lg_n_',n_int,'_ir_',inv_r,'_pX_',p_X,'.rds'))
                }
            }
        }}}
