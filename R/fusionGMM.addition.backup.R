library(glmnet)
pseudo_Xy_multiv_addition<-function(
        C_half,X,A,y,study_info){
    XA<-cbind(X,A)
    xatxa<-crossprod(XA)
    xtxa<-crossprod(X,XA)
    #pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_X<-C_half%*%rbind(xatxa,xtxa)
    #pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    pseudo_y1<-crossprod(y,XA)
    #pseudo_y2<-crossprod(UKBB_pop[,var_SNP],UKBB_pop[,var_SNP])%*%study_info[[1]]$Coeff
    pseudo_y2<-crossprod(X)%*%study_info[[1]]$Coeff
    pseudo_y<-c(c(pseudo_y1,pseudo_y2)%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}
# work for computation of weighting matrix C
var_U_beta_theta_func_multiv_addition<-function(
        X,A,y,beta,
        study_info){
    XA<-cbind(X,A)
    #xtbeta<-(UKBB_pop[,-1]%*%beta)
    xatbeta<-(XA%*%beta)
    #xttheta<-(UKBB_pop[,var_SNP]%*%study_info[[1]]$Coeff)
    xttheta<-(X%*%study_info[[1]]$Coeff)
    #var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    var_11<-crossprod(XA*c(xatbeta-y))
    #var_22<-crossprod(UKBB_pop[,var_SNP]*c(xtbeta-xttheta),UKBB_pop[,var_SNP]*c(xtbeta-xttheta))
    var_22<-crossprod(X*c(xatbeta-xttheta))
    #var_12<-crossprod(UKBB_pop[,-1]*c(xtbeta-UKBB_pop[,1]),UKBB_pop[,var_SNP]*c(xtbeta-xttheta))
    var_12<-crossprod(XA*c(xatbeta-y),X*c(xatbeta-xttheta))
    (1/nrow(X))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}

pseudo_Xy_univ_addition<-function(
        C_half,X,A,y,study_info){
    XA<-cbind(X,A)
    xatxa<-crossprod(XA)
    xtxa<-crossprod(X,XA)
    #pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_X<-C_half%*%rbind(xatxa,xtxa)
    #pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    pseudo_y1<-crossprod(y,XA)
    #pseudo_y22<-sapply(1:N_SNP, function(snp_id){
    #    u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
    #    c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    #})
    pseudo_y2<-sapply(1:ncol(X), function(id){
        u2_id<-c(X[,id]*study_info[[id]]$Coeff)
        c(u2_id%*%X[,id])
    })
    pseudo_y<-c(c(pseudo_y1,pseudo_y2)%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}
# work for computation of weighting matrix C
var_U_beta_theta_func_univ_addition<-function(
        X,A,y,beta,
        study_info){
    XA<-cbind(X,A)
    #xtbeta<-(UKBB_pop[,-1]%*%beta)
    xatbeta<-(XA%*%beta)
    #var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    var_11<-crossprod(XA*c(xatbeta-y))
    #u2_theta_coef<-sapply(1:ncol(X), function(id){
    #    #xtgamma_id<-c(UKBB_pop[,paste0("SNP",snp_id)]*study_info[[snp_id]]$Coeff)
    #    xtgamma_id<-c(X[,id]*study_info[[id]]$Coeff)
    #    xtbeta-xtgamma_id
    #}) #col is SNP #row is sample

    u2_theta_coef<-sapply(1:ncol(X), function(id){
        xtgamma_id<-c(X[,id]*study_info[[id]]$Coeff)
        xatbeta-xtgamma_id
    }) #row is sample #col is SNP
    #var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
    var_22<-crossprod(u2_theta_coef*X)
    #var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),u2_theta_coef*UKBB_pop[,var_SNP])
    var_12<-crossprod(XA*c(xatbeta-y),u2_theta_coef*X)
    (1/nrow(X))*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}


fusionGMM.addition<-function(
        X,A,y,
        study_info,
        summary_type = "multi",
        penalty_type = "lasso",
        initial_with_type = "ridge",
        initial_with_GMM = FALSE,
        remove_penalty_X = FALSE,
        remove_penalty_A = FALSE,
        tune_ratio = TRUE,
        fix_ratio = NULL,
        ratio_lower = NULL,
        ratio_upper = NULL,
        ratio_count = 10,
        gamma_adaptivelasso = 1/2,
        inference = FALSE,
        approx_cross_validation = TRUE,
        kfolds = 10,
        desparseC = FALSE
){
    if (!summary_type %in% c("uni","multi")){
        stop("Select summary_type for input summary statistics(study_info).
       Use 'uni' for univariate summary statistics.
       Use 'multi' for multivariate summary statistics.")
    }
    if(!penalty_type %in% c("adaptivelasso","lasso","ridge")){
        stop("Select penalty type from c('adaptivelasso','lasso','ridge').")
    }
    final_alpha = 1
    if(penalty_type == "ridge"){final_alpha = 0}
    if (summary_type == "uni"){
        pseudo_Xy <- pseudo_Xy_univ_addition
        var_U_beta_theta_func <- var_U_beta_theta_func_univ_addition
    }else if (summary_type == "multi"){
        pseudo_Xy <- pseudo_Xy_multiv_addition
        var_U_beta_theta_func <- var_U_beta_theta_func_multiv_addition
    }

    if(!is.null(fix_ratio)){
        if(tune_ratio){
            stop("If ratio is fixed, please set tune_ratio as FALSE")
        }else if(remove_penalty_X | remove_penalty_A){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }
    # X and y should be scaled to avoid potential mismatch of intercept term.
    # n for sample size; p for feature size
    nX<-nrow(X)
    nXext<-study_info[[1]]$Sample_size
    pX<-ncol(X)
    pA<-ncol(A)
    if(nX<pX+pA){desparseC<-TRUE}
    if(nX<700){desparseC<-TRUE}
    if(nX<700){initial_with_GMM<-TRUE}
    XA<-cbind(X,A)
    Xid<-1:pX
    Aid<-(pX+1):(pX+pA)
    fix_penalty<-rep(1,pX+pA)
    if(remove_penalty_X){
        fix_penalty[Xid]<-0
    }
    if(remove_penalty_A){
        fix_penalty[Aid]<-0
    }
    if(!is.null(fix_ratio)){
        fix_penalty[Xid]<-fix_ratio
    }
    xatxa<-crossprod(XA)
    xtx<-crossprod(X)


    if (summary_type == "uni"){
        diag_theta_sd<-sapply(1:pX,function(id){
            sqrt(study_info[[id]]$Covariance)*c(X[,id]%*%X[,id])})
        var_U2<-t(diag_theta_sd*cor(X))*diag_theta_sd/nX
        #C_22_half<-diag(sapply(1:pX,function(id){
        #    sqrt(study_info[[1]]$Sample_size/(study_info[[id]]$Covariance))/c(X[,id]%*%X[,id])
        #}))
    }else if (summary_type == "multi"){
        var_U2<-xtx%*%study_info[[1]]$Covariance%*%xtx/nX
    }
    if(initial_with_type %in% c("ridge","lasso")){
        if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
        fit_initial<-cv.glmnet(x=XA,y= y,alpha = initial_alpha,penalty.factor = fix_penalty)
        beta_initial<-c(coef(fit_initial,s="lambda.min")[-1])
        if (penalty_type == "adaptivelasso"){
            w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
            w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
            w_adaptive<-w_adaptive*fix_penalty
        }else{
            w_adaptive<-fix_penalty
        }
    }else{
        # Initial estimation of C
        var_U1<-xatxa/nX*c(var(y))
        if(desparseC){
            C_11_half<-diag(1/sqrt(diag(var_U1)))
            C_22_half<-diag(1/sqrt(diag(var_U2)))
        }else{
            inv_C_11_svd<-corpcor::fast.svd(var_U1)
            C_11_half<-inv_C_11_svd$u%*%diag(1/sqrt(inv_C_11_svd$d))%*%t(inv_C_11_svd$v)
            inv_C_22_svd<-corpcor::fast.svd(var_U2)
            C_22_half<-inv_C_22_svd$u%*%diag(1/sqrt(inv_C_22_svd$d))%*%t(inv_C_22_svd$v)
        }

        C_half<-adiag(C_11_half,C_22_half)

        # Prepare for initial model
        pseudo_Xy_list<-pseudo_Xy(C_half,X,A,y,study_info)
        initial_sf<-nX/sqrt(nrow(pseudo_Xy_list$pseudo_X))
        pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
        pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
        fit_initial<-cv.glmnet(x= pseudo_X,y= pseudo_y,standardize=F,intercept=F,alpha = 0,penalty.factor = fix_penalty)
        beta_initial<-c(coef(fit_initial,s ='lambda.min')[-1])
        if (penalty_type == "adaptivelasso"){
            #fit_for_weight<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,intercept=F,lambda = 0.01,alpha = 0.01)
            #beta_for_weight<-c(coef(fit_for_weight)[-1])
            #w_adaptive<-1/abs(beta_for_weight)^gamma_adaptivelasso
            #w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
            #w_adaptive<-w_adaptive*fix_penalty
            w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
            w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
            w_adaptive<-w_adaptive*fix_penalty
        }else{
            w_adaptive<-fix_penalty
        }
    }

    # Refined estimation of C
    var_1st_U_beta_theta<-var_U_beta_theta_func(X = X,A = A, y = y,beta = beta_initial,
                                                                study_info = study_info)
    var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = pX+pA,ncol = pX+pA),var_U2)
    inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat
    if(desparseC){
        C_half<-diag(1/sqrt(diag(inv_C)))
    }else{
        inv_C_svd<-corpcor::fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
        C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
    }

    # Prepare for final model
    pseudo_Xy_list<-pseudo_Xy(C_half,X,A,y,study_info)
    initial_sf<-nX/sqrt(nrow(pseudo_Xy_list$pseudo_X))
    pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
    pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf

    # Fit final model
    fit_final<-cv.glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                         intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
    lambda_list<-fit_final$lambda

    if(tune_ratio & !remove_penalty_X & !remove_penalty_A){
        if(is.null(ratio_lower)){ratio_lower<-sqrt(nX/(nX+nXext))/2}
        if(is.null(ratio_upper)){ratio_upper<-1}
        ratio_range<-exp(seq(log(ratio_lower),log(ratio_upper),(log(ratio_upper)-log(ratio_lower))/ratio_count))
        ratio_range<-c(ratio_range,2)
        if(approx_cross_validation){
            ## approximate cross validation
            cv_ratio<-lapply(ratio_range,function(cur_ratio){
                ratio_vec<-c(rep(cur_ratio,pX),rep(1,pA))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv.glmnet(x=pseudo_X,y=pseudo_y,
                          standardize=F,intercept=F,alpha = final_alpha,
                          penalty.factor=w_adaptive_ratio,
                          lambda = lambda_list)
            })
            cv_ratio_mse<-sapply(1:length(ratio_range), function(i){
                min(cv_ratio[[i]]$cvm)})
            fit_final_ratio<-cv_ratio[[which.min(cv_ratio_mse)]]
            beta<-c(coef(fit_final_ratio,s="lambda.min")[-1])
            return_list<-list("beta"=beta,
                              "lambda_list"=lambda_list,
                              "ratio_list"=ratio_range,
                              "lambda_min"=fit_final_ratio$lambda.min,
                              "ratio_min"=ratio_range[which.min(cv_ratio_mse)])
        }else{
            ## detailed cross validation
            if(length(unique(y)) <= 2){
                index_fold<-caret::createFolds(as.numeric(y>0),k = kfolds)
            }else{
                index_fold<-caret::createFolds(y,k = kfolds)
            }

            mse_lam_ratio<-lapply(1:kfolds, function(cur_fold){
                index_test<-index_fold[[cur_fold]]
                Xtrain<-X[-index_test,]
                Xtest<-X[index_test,]
                Atrain<-A[-index_test,]
                Atest<-A[index_test,]
                ytrain<-y[-index_test]
                ytest<-y[index_test]
                pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,
                                                                ytrain,study_info)
                initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
                pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
                pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
                mse_lam_ratio_fold<-sapply(lambda_list,function(cur_lam){
                    sapply(ratio_range,function(cur_ratio){
                        ratio_vec<-c(rep(cur_ratio,pX),rep(1,pA))
                        w_adaptive_ratio<-w_adaptive*ratio_vec
                        cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                                       standardize=F,intercept=F,alpha = final_alpha,
                                       penalty.factor = w_adaptive_ratio,
                                       lambda = cur_lam)
                        cur_beta<-coef(cv_fit)[-1]
                        mean((cbind(Xtest,Atest)%*%cur_beta - ytest)^2)
                    })
                }) # row is ratio_range & col is lambda_list
                mse_lam_ratio_fold
            })

            sum_mse_lam_ratio<-Reduce(`+`, mse_lam_ratio)
            ids<-which(sum_mse_lam_ratio==min(sum_mse_lam_ratio),arr.ind = TRUE)
            final.ratio.min<-ratio_range[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
            ratio_vec<-c(rep(final.ratio.min,pX),rep(1,pA))
            w_adaptive_ratio<-w_adaptive*ratio_vec
            fit_final_lam_ratio<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                        intercept=F,alpha = final_alpha,
                                        penalty.factor = w_adaptive_ratio,
                                        lambda = final.lambda.min)

            beta<-coef(fit_final_lam_ratio)[-1]
            return_list<-list("beta"=beta,
                              "lambda_list"=lambda_list,
                              "ratio_list"=ratio_range,
                              "mse_fold"=mse_lam_ratio,
                              "lambda_min"=final.lambda.min,
                              "ratio_min"=final.ratio.min)
        }

    }else{
        # NOT tune ratio
        if(approx_cross_validation){
            # approximate cross validation
            beta<-coef(fit_final,s ='lambda.min')[-1]
            return_list<-list("beta"=beta,
                              "lambda_list"=lambda_list,
                              "lambda_min"=fit_final$lambda.min)
        }else{
            ## detailed cross validation
            if(length(unique(y)) <= 2){
                index_fold<-caret::createFolds(as.numeric(y>0),k = kfolds)
            }else{
                index_fold<-caret::createFolds(y,k = kfolds)
            }
            mse_fold<-sapply(1:kfolds, function(cur_fold){
                index_test<-index_fold[[cur_fold]]
                Xtrain<-X[-index_test,]
                Xtest<-X[index_test,]
                Atrain<-A[-index_test,]
                Atest<-A[index_test,]
                ytrain<-y[-index_test]
                ytest<-y[index_test]
                pseudo_Xy_list_train<-pseudo_Xy(C_half,Xtrain,Atrain,
                                                                ytrain,study_info)
                initial_sf_train<-nrow(Xtrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
                pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
                pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
                mse_lam<-sapply(lambda_list,function(cur_lam){
                    cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                                   standardize=F,intercept=F,
                                   alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
                    cur_beta<-coef(cv_fit)[-1]
                    sum((cbind(Xtest,Atest)%*%cur_beta - ytest)^2)
                })
                mse_lam
            })
            final.lambda.min<-lambda_list[which.min(rowSums(mse_fold))]
            fit_final_cv<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                 intercept=F,alpha = final_alpha,
                                 penalty.factor = w_adaptive,
                                 lambda = final.lambda.min)
            beta<-c(coef(fit_final_cv,s ='lambda.min')[-1])
            return_list<-list("beta"=beta,
                              "lambda_list"=lambda_list,
                              "mse_fold"=mse_fold,
                              "lambda_min"=final.lambda.min)
        }
    }

    index_nonzero<-which(beta!=0)
    if(inference & length(index_nonzero) > 1){
        Sigsum_half<-cbind(xatxa/nX,crossprod(XA,X)/nX)%*%C_half
        Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
        Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
        inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
        final_v<-diag(inv_Sigsum_scaled_nonzero)

        pval_final<-pchisq(nX*beta[index_nonzero]^2/final_v,1,lower.tail = F)
        corrected_pos<-index_nonzero[which(pval_final<0.05/length(index_nonzero))]
        return_list<-c(return_list,
                       list("corrected_pos"=pos,
                            "nonzero_pos"=index_nonzero,
                            "pval"=pval_final,
                            "nonzero_var"=final_v))
    }
    return(return_list)
}
