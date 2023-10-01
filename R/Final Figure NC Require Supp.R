library(ggplot2)
library(latex2exp)

theme_bar <- function(..., bg='white'){
    require(grid)
    theme_classic(...) +
        theme(rect=element_rect(fill=bg),
              plot.margin=unit(rep(0.5,4), 'lines'),
              panel.background=element_rect(fill='transparent', color='black'),
              panel.border=element_rect(fill='transparent', color='transparent'),
              panel.grid=element_blank(),
              axis.title.y=element_text(face = "bold",size = 14),
              axis.text = element_text(size = 9),
              axis.text.x = element_text(face = "bold"),
              axis.text.y = element_text(face = "bold"),
              legend.text = element_text(face = "bold",size = 9),
              legend.background = element_rect( linetype="solid",colour ="black"),
              plot.title = element_text(face = "bold")
        )

}
p_X_list<-c(10,40)
inv_r_list<-c(10,30)
#n_int_list<-floor(exp(seq(log(100),log(10000),(log(10000) - log(100))/9)))
n_int_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)[2:10]
n_int_list2<-n_int_list
n_int_list2[1]<-220
lnlg<-"Lg"
lnlg2<-"Lg"#"LgAP"
AUCR2<-"AUC"
gname="Lasso"
gname="AdaptiveLasso"
extra_fold<-"_same"
filecount=paste0("_",2)
#filecount=""
####extra_fold<-""
####extra_fold<-"_sC2"
lowbound1<-1
lowbound2<-1
upbound1<-0
upbound2<-0
for (id1 in 1:2){
    for(id2 in 1:2){
        p_X<-p_X_list[id1]
        inv_r<-inv_r_list[id2]
        rrid3<-c()
        bgid3<-c()
        skipid<-c()
        message(c(id1,id2))
        for(id3 in c(1:length(n_int_list))){
            p_X<-p_X_list[id1]
            inv_r = inv_r_list[id2]
            n_int = n_int_list[id3]
            rr<-c()
            bg<-c()
            if(gname == "Lasso"){
                lassoindex<-c(4,5,9,10,13,14)
                lassoindex1<-c(1,5,6,9,10)
            }else{
                lassoindex<-c(4,7,11,12,15,16)
                lassoindex1<-c(3,7,8,11,12)
            }

            res<-readRDS(paste0('/users/rzhao1/intgmm/result',extra_fold,'/',lnlg,'/',lnlg2,'_n_',n_int,'_ir_',inv_r,'_pX_',p_X,filecount,'.rds'))
            for(i in 1:length(res)){
                rr<-rbind(rr,res[[i]]$rr)
            }
            highbase<-0.754
            #highbase<-0.343
            rr<-rr[,lassoindex]
            if(AUCR2 == "AUC"){
                if(lnlg == "Lg" & id3 == 1){
                    #if(id1==1){
                    #    rm_index<-c(8,14)
                    #}else{
                    #    rm_index<-c(94)
                    #}
                    #print(rr[rm_index,])
                    #print(bg[rm_index,])
                    #print(which(rowSums(bg[,-1]>3)>0))
                    #rr<-rr[rowSums(bg[,-1]>3)==0,]
                    #bg<-bg[rowSums(bg[,-1]>3)==0,]
                    #rr<-rr[-rm_index,]
                    #bg<-bg[-rm_index,]
                }
            }

            rrid3<-rbind(rrid3,colMeans(rr))

        }


        rrid3 = rrid3[,c(1,2,3,4)]

        cm<-function(methods){
            methods2<-methods
            methods2[which(methods== "rmul")]<-"HTLGMM w/ Lasso2"
            methods2[which(methods== "rmul1")]<-"HTLGMM w/ Lasso"
            methods2[which(methods== "rmulada")]<-"HTLGMM w/ Lasso2"
            methods2[which(methods== "rmulada1")]<-"HTLGMM w/ Lasso"
            methods2[which(substr(methods,1,4)%in%c("rlas"))]<-"Internal Only w/ Lasso"
            methods2[which(substr(methods,1,4)%in%c("rada"))]<-"Internal Only w/ Lasso"
            methods2[which(substr(methods,1,2)%in%c("rX"))]<-"ExternalOnly"
            list("m1"=methods,"m2"=methods2)
        }
        cmrr<-cm(colnames(rrid3))

        #print(colnames(rrid3))
        lowbound1<-min(lowbound1,rrid3)

        upbound1<-max(upbound1,rrid3)
        #print(rrid3)

        df<-data.frame(SqrtSize=rep(sqrt(n_int_list),ncol(rrid3)),
                       Method=rep(cmrr$m2,1,each=nrow(rrid3)),
                       R2=c(matrix(rrid3,nrow=1,byrow = T)))
        if(gname == "Lasso"){
            df$Method<-factor(df$Method, levels = c("HTLGMM w/ Lasso","Internal Only w/ Lasso","ExternalOnly","HTLGMM w/ Lasso2"))
        }else{
            df$Method<-factor(df$Method, levels = c("HTLGMM w/ Lasso","Internal Only w/ Lasso","ExternalOnly","HTLGMM w/ Lasso2"))
        }

        #print(df)
        library(ggplot2)
        #n_int_list2[1]<-400
        n_x<-c(200,400,600,900,1200,1600,2000,2500,3000)
        gg<-ggplot(data = df)+
            geom_point(aes(x=SqrtSize,y=R2,color=Method,shape = Method))+
            geom_line(aes(x=SqrtSize,y=R2,color=Method,linetype = Method))+
            geom_hline(yintercept = highbase,color = "darkred")+
            #labs(title = paste0(AUCR2,":Comp:Incomp=1:",inv_r,";RiskFactor:",p_X))+
            #ylab(AUCR2)+
            theme_bar()+
            scale_x_continuous(breaks = sqrt(n_int_list),labels =n_x,name = "Internal Sample Size")

        if(AUCR2=="AUC"){
            if(length(n_int_list) == 9){
                gg<-gg+scale_y_continuous(limits = c(0.66, 0.755))
            }

            gg<-gg+annotate("text", x = 33, y = highbase*1.005, label = "",size=3)+ylab(TeX(r'($AUC$)'))
        }else{

            ## R2 R2 R2 R2
            if(length(n_int_list) == 9){
                #gg<-gg+scale_y_continuous(limits = c(0.205, 0.355))
                gg<-gg+scale_y_continuous(limits = c(0.18, 0.355))
            }
            gg<-gg+annotate("text", x = 33, y = highbase*1.015, label = "",size=3)+ylab(TeX(r'($R^2$)'))
        }
        if(p_X==10 & inv_r == 10){
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic AUC: $p_X=10,p_A=150,n_{ext}/n=10$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear $R^2$: $p_X=10,p_A=150,n_{ext}/n=10$)'))
            }
        }else if(p_X==10 & inv_r == 30){
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic AUC: $p_X=10,p_A=150,n_{ext}/n=30$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear $R^2$: $p_X=10,p_A=150,n_{ext}/n=30$)'))
            }
        }else if(p_X==40 & inv_r == 10){
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic AUC: $p_X=40,p_A=150,n_{ext}/n=10$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear $R^2$: $p_X=40,p_A=150,n_{ext}/n=10$)'))
            }
        }else{
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic AUC: $p_X=40,p_A=150,n_{ext}/n=30$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear $R^2$: $p_X=40,p_A=150,n_{ext}/n=30$)'))
            }
        }
        gg<-gg+scale_color_manual(values = c("#F8766D","#00BA38","#619CFF","#F564E3"))
        # Use this palette in your plot
        #   gg<-gg+scale_color_manual(values = c("#F8766D","#00BA38","#619CFF","#F564E3"))
        #ggsave(paste0("/users/rzhao1/intgmm/plot/plot",lnlg,extra_fold,"/",lnlg2,gname,":",AUCR2,":Comp:Incomp=1:",inv_r,";RiskFactor:",p_X,".png"),gg,width = 7,height = 4.5)
        ggsave(paste0("/users/rzhao1/intgmm/finalplot/dbplot",lnlg2,gname,":",AUCR2,":Comp:Incomp=1:",inv_r,";RiskFactor:",p_X,".png"),gg,width = 7,height = 4.5)
    }
}

