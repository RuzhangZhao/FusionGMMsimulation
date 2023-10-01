library(ggplot2)
library(latex2exp)

theme_bar <- function(..., bg='white'){
    require(grid)
    theme_classic(...) +
        theme(rect=element_rect(fill=bg),
              plot.margin=unit(rep(0.5,4), 'lines'),
              panel.background=element_rect(fill='transparent', color='black'),
              panel.border=element_rect(fill='transparent', color='transparent'),
              panel.grid=element_blank(),#去网格线
              #axis.title.x = element_blank(),#去x轴标签
              axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
              axis.text = element_text(size = 9),#坐标轴刻度标签加粗
              axis.text.x = element_text(face = "bold"),
              axis.text.y = element_text(face = "bold"),
              # axis.ticks = element_line(color='black'),#坐标轴刻度线
              # axis.ticks.margin = unit(0.8,"lines"),
              #legend.title=element_blank(),#去除图例标题
              #legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
              #legend.position=c(0.5, 0.7),#图例在绘图区域的位置
              #legend.position='top',#图例放在顶部
              #legend.direction = "horizontal",#设置图例水平放置
              # legend.spacing.x = unit(2, 'cm'),
              legend.text = element_text(face = "bold",size = 9),
              legend.background = element_rect( linetype="solid",colour ="black"),
              plot.title = element_text(face = "bold")
              # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
              # legend.box.margin =margin(-10,0,0,0)
        )

}
p_X_list<-c(10,40)
inv_r_list<-c(10,30)
#n_int_list<-floor(exp(seq(log(100),log(10000),(log(10000) - log(100))/9)))
n_int_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)[2:10]
n_int_list2<-n_int_list
n_int_list2[1]<-220
lnlg<-"Ln"
lnlg2<-"Ln"#"LgAP"
AUCR2<-"R2"
gname="Lasso"
#gname="AdaptiveLasso"
extra_fold<-"_same"
filecount=paste0("_",2)
#filecount=""
#extra_fold<-""
#extra_fold<-"_sC2"
lowbound1<-1
lowbound2<-1
upbound1<-0
upbound2<-0
for (id1 in c(1,2)){
    p_X<-p_X_list[id1]
    for(id2 in c(1,2)){
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

            if(!file.exists(paste0('/users/rzhao1/intgmm/result',extra_fold,'/',lnlg,'/',lnlg2,'_n_',n_int,'_ir_',inv_r,'_pX_',p_X,filecount,'.rds'))){
                skipid<-c(skipid,id3)
            }else{
                res<-readRDS(paste0('/users/rzhao1/intgmm/result',extra_fold,'/',lnlg,'/',lnlg2,'_n_',n_int,'_ir_',inv_r,'_pX_',p_X,filecount,'.rds'))
                for(i in 1:length(res)){
                    rr<-rbind(rr,res[[i]]$rr)
                    bg<-rbind(bg,res[[i]]$beta)
                }
                highbase<-0.754
                bg<-bg[,lassoindex1]
                rr<-rr[,lassoindex]
                if(AUCR2 == "AUC"){
                    if(lnlg == "Lg" & id3 == 1){
                        if(id1==1){
                            rm_index<-c(8,14)
                        }else{
                            rm_index<-c(94)
                        }
                        print(rr[rm_index,])
                        print(bg[rm_index,])
                        #print(which(rowSums(bg[,-1]>3)>0))
                        #rr<-rr[rowSums(bg[,-1]>3)==0,]
                        #bg<-bg[rowSums(bg[,-1]>3)==0,]
                        rr<-rr[-rm_index,]
                        bg<-bg[-rm_index,]
                    }
                }
                #if(id3<2){print(bg)}

                rrid3<-rbind(rrid3,colMeans(rr))
                bgid3<-rbind(bgid3,colMeans(bg))

            }}
        #print(bgid3)
        #[1] "r0"       "r0X"      "r0A"      "rX"       "rlasso"   "rlasso1"  "rada"     "rada1"
        #[9] "rmul"     "rmul1"    "rmulada"  "rmulada1" "runi"     "runi1"    "runiada"  "runiada1"
        cm<-function(methods){
            methods2<-methods
            methods2[which(substr(methods,1,4)%in%c("rmul"))]<-"TLGMM-multivariate"
            methods2[which(substr(methods,1,4)%in%c("runi"))]<-"TLGMM-univariate"
            methods2[which(substr(methods,1,4)%in%c("rlas","rada"))]<-"InternalOnly"
            methods2[which(substr(methods,1,2)%in%c("rX"))]<-"ExternalOnly"
            methods[which(methods%in%c("rmul1","runi1"))]<-"TLGMM-Lasso-1lam"
            methods[which(methods%in%c("rmul","runi"))]<-"TLGMM-Lasso-2lam"
            methods[which(methods%in%c("rmulada1","runiada1"))]<-"TLGMM-AdaLasso-1lam"
            methods[which(methods%in%c("rmulada","runiada"))]<-"TLGMM-AdaLasso-2lam"
            methods[which(methods%in%c("rlasso"))]<-"naive-Lasso"
            methods[which(methods%in%c("rlasso1"))]<-"naive-Lasso-lambda.1se"
            methods[which(methods%in%c("rada"))]<-"naive-AdaLasso"
            methods[which(methods%in%c("rada1"))]<-"naive-AdaLas-lambda.1se"
            methods[which(methods%in%c("rX"))]<-"OLS"
            list("m1"=methods,"m2"=methods2)
        }
        cmrr<-cm(colnames(rrid3))
        cmbg<-cm(colnames(bgid3))
        #print(colnames(rrid3))
        lowbound1<-min(lowbound1,rrid3)
        lowbound2<-min(lowbound2,bgid3)
        upbound1<-max(upbound1,rrid3)
        upbound2<-max(upbound2,bgid3)
        #print(rrid3)
        if(!is.null(skipid)){
            print("skip")
            df<-data.frame(SqrtSize=rep(sqrt(n_int_list)[-skipid],ncol(rrid3)),
                           Method=rep(cmrr$m1,1,each=nrow(rrid3)),
                           DataSource=rep(cmrr$m2,1,each=nrow(rrid3)),
                           R2=c(matrix(rrid3,nrow=1,byrow = T)))
            df2<-data.frame(SqrtSize=rep(sqrt(n_int_list)[-skipid],ncol(bgid3)),
                            Method=rep(cmbg$m1,1,each=nrow(bgid3)),
                            DataSource=rep(cmbg$m2,1,each=nrow(bgid3)),
                            BetaMSE=c(matrix(bgid3,nrow=1,byrow = T)))
            df$DataSource<-factor(df$DataSource, levels = c("TLGMM-multivariate","TLGMM-univariate","InternalOnly"))
            df2$DataSource<-factor(df2$DataSource, levels = c("TLGMM-multivariate","TLGMM-univariate","InternalOnly"))
        }else{
            df<-data.frame(SqrtSize=rep(sqrt(n_int_list),ncol(rrid3)),
                           Method=rep(cmrr$m1,1,each=nrow(rrid3)),
                           DataSource=rep(cmrr$m2,1,each=nrow(rrid3)),
                           R2=c(matrix(rrid3,nrow=1,byrow = T)))
            df2<-data.frame(SqrtSize=rep(sqrt(n_int_list),ncol(bgid3)),
                            Method=rep(cmbg$m1,1,each=nrow(bgid3)),
                            DataSource=rep(cmbg$m2,1,each=nrow(bgid3)),
                            BetaMSE=c(matrix(bgid3,nrow=1,byrow = T)))
            if(gname == "Lasso"){
                df$Method<-factor(df$Method, levels = c("TLGMM-Lasso-1lam","TLGMM-Lasso-2lam","naive-Lasso","OLS"))
                df2$Method<-factor(df2$Method, levels = c("TLGMM-Lasso-1lam","TLGMM-Lasso-2lam","naive-Lasso","OLS"))
            }else{
                df$Method<-factor(df$Method, levels = c("TLGMM-AdaLasso-1lam","TLGMM-AdaLasso-2lam","naive-AdaLasso","OLS"))
                df2$Method<-factor(df2$Method, levels = c("TLGMM-AdaLasso-1lam","TLGMM-AdaLasso-2lam","naive-AdaLasso","OLS"))
            }

            df$DataSource<-factor(df$DataSource, levels = c("TLGMM-multivariate","TLGMM-univariate","InternalOnly","ExternalOnly"))

            df2$DataSource<-factor(df2$DataSource, levels = c("TLGMM-multivariate","TLGMM-univariate","InternalOnly"))
        }
        #print(df)
        library(ggplot2)
        #n_int_list2[1]<-400
        gg<-ggplot(data = df)+
            geom_point(aes(x=SqrtSize,y=R2,color=Method,shape = DataSource))+
            geom_line(aes(x=SqrtSize,y=R2,color=Method,linetype = DataSource))+
            geom_hline(yintercept = highbase,color = "darkred")+
            #labs(title = paste0(AUCR2,":Comp:Incomp=1:",inv_r,";RiskFactor:",p_X))+
            #ylab(AUCR2)+
            theme_bar()+
            scale_x_continuous(breaks = sqrt(n_int_list),labels =n_int_list2,name = "Internal Sample Size")

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

        my_palette <- c("TLGMM-multivariate" = "#F8766D",
                        "TLGMM-univariate" = "#00BA38",
                        "InternalOnly" = "#619CFF",
                        "ExternalOnly" = "#F564E3")

        # Use this palette in your plot
        gg<-gg+scale_color_manual(values = c("#F8766D","#00BA38","#619CFF","#F564E3"))
        ggsave(paste0("/users/rzhao1/intgmm/plot2/plot",lnlg,extra_fold,"/",lnlg2,gname,":",AUCR2,":Comp:Incomp=1:",inv_r,";RiskFactor:",p_X,".png"),gg,width = 7,height = 4.5)
        gg<-ggplot(data = df2)+
            geom_point(aes(x=SqrtSize,y=BetaMSE,color=Method,shape = DataSource))+
            geom_line(aes(x=SqrtSize,y=BetaMSE,color=Method,linetype = DataSource))+
            theme_bar()+
            scale_x_continuous(breaks = sqrt(n_int_list),labels =n_int_list2,name = "Internal Sample Size")+
            ylab(TeX(r'(Estimation Error of $\beta$:)'))
        if(AUCR2 == "AUC"){
            if(length(n_int_list) == 9){
                #gg<-gg+scale_y_continuous(limits = c(0.08, 0.52))
                gg<-gg+scale_y_continuous(limits = c(0.09, 0.937))
            }

        }else{
            if(length(n_int_list) == 9){
                #gg<-gg+scale_y_continuous(limits = c(0.15, 1.64))
                gg<-gg+scale_y_continuous(limits = c(0.155, 2.33))
            }
        }
        if(p_X==10 & inv_r == 10){
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic: $p_X=10,p_A=150,n_{ext}/n=10$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear: $p_X=10,p_A=150,n_{ext}/n=10$)'))
            }
        }else if(p_X==10 & inv_r == 30){
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic: $p_X=10,p_A=150,n_{ext}/n=30$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear: $p_X=10,p_A=150,n_{ext}/n=30$)'))
            }
        }else if(p_X==40 & inv_r == 10){
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic: $p_X=40,p_A=150,n_{ext}/n=10$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear: $p_X=40,p_A=150,n_{ext}/n=10$)'))
            }
        }else{
            if(AUCR2=="AUC"){
                gg<-gg+labs(title = TeX(r'(Logistic: $p_X=40,p_A=150,n_{ext}/n=30$)'))
            }else{
                gg<-gg+labs(title = TeX(r'(Linear: $p_X=40,p_A=150,n_{ext}/n=30$)'))
            }
        }
        ggsave(paste0("/users/rzhao1/intgmm/plot2/plot",lnlg,extra_fold,"/",lnlg2,gname,":BetaMSE:Comp:Incomp=1:",inv_r,";RiskFactor:",p_X,".png"),gg,width = 7,height = 4.5)
    }
}
print(c(lowbound1,upbound1,lowbound2,upbound2))

