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
              axis.text = element_text(face = "bold",size = 9),#坐标轴刻度标签加粗
              # axis.ticks = element_line(color='black'),#坐标轴刻度线
              # axis.ticks.margin = unit(0.8,"lines"),
              #legend.title=element_blank(),#去除图例标题
              #legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
              #legend.position=c(0.5, 0.7),#图例在绘图区域的位置
              #legend.position='top',#图例放在顶部
              #legend.direction = "horizontal",#设置图例水平放置
              # legend.spacing.x = unit(2, 'cm'),
              legend.text = element_text(face = "bold",size = 9),
              legend.background = element_rect( linetype="solid",colour ="black")
              # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
              # legend.box.margin =margin(-10,0,0,0)
        )

}

p_X_list<-c(10,40)
inv_r_list<-c(10,30)
#n_int_list<-floor(exp(seq(log(100),log(10000),(log(10000) - log(100))/9)))
n_int_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)
n_int_list
for (id1 in c(2)){
    p_X<-p_X_list[id1]
    for(id2 in c(2)){
        inv_r<-inv_r_list[id2]
        rrid3<-c()
        bgid3<-c()
        skipid<-c()
        message(c(id1,id2))
        for(id3 in c(1:6)){
            p_X<-p_X_list[id1]
            inv_r = inv_r_list[id2]
            n_int = n_int_list[id3]
            rr<-c()
            bg<-c()
            gname="Lasso"
            #gname="AdaptiveLasso"
            if(gname == "Lasso"){
                lassoindex<-c(5,6,9,10,13,14)
                lassoindex1<-c(1,2,5,6,9,10)
            }else{
                lassoindex<-c(7,8,11,12,15,16)
                lassoindex1<-c(3,4,7,8,11,12)
            }


            if(!file.exists(paste0('../FusionGMMdata/Ln/Ln_n_',n_int,'_ir_',inv_r,'_pX_',p_X,'.rds'))){
                skipid<-c(skipid,id3)
            }else{
                res<-readRDS(paste0('../FusionGMMdata/Ln/Ln_n_',n_int,'_ir_',inv_r,'_pX_',p_X,'.rds'))
                for(i in 1:length(res)){
                    rr<-rbind(rr,res[[i]]$rr)
                    bg<-rbind(bg,res[[i]]$beta)
                }
                bg<-bg[,lassoindex1]
                rr<-rr[,lassoindex]
                bg<-bg[rowSums(rr<0)==0,]
                rr<-rr[rowSums(rr<0)==0,]
                #if(id3%in%c(2)){rr<-rr[-6,]
                 #   bg<-bg[-6,]}
                rrid3<-rbind(rrid3,colMeans(rr))
                bgid3<-rbind(bgid3,colMeans(bg))
            }
        }
        highbase<-rrid3[1,1]
        #[1] "r0"       "r0X"      "r0A"      "rX"       "rlasso"   "rlasso1"  "rada"     "rada1"
        #[9] "rmul"     "rmul1"    "rmulada"  "rmulada1" "runi"     "runi1"    "runiada"  "runiada1"
        cm<-function(methods){
            methods2<-methods
            methods2[which(substr(methods,1,4)%in%c("rmul"))]<-"IntGMM-multivariate"
            methods2[which(substr(methods,1,4)%in%c("runi"))]<-"IntGMM-univariate"
            methods2[which(substr(methods,1,4)%in%c("rlas","rada"))]<-"InternalOnly"
            methods2[which(substr(methods,1,2)%in%c("rX"))]<-"ExternalOnly"
            methods[which(methods%in%c("rmul1","runi1"))]<-"IntGMM-Lasso-1lam"
            methods[which(methods%in%c("rmul","runi"))]<-"IntGMM-Lasso-2lam"
            methods[which(methods%in%c("rmulada1","runiada1"))]<-"IntGMM-AdaLasso-1lambda"
            methods[which(methods%in%c("rmulada","runiada"))]<-"IntGMM-AdaLasso-2lambda"
            methods[which(methods%in%c("rlasso"))]<-"naive-Lasso-lambda.min"
            methods[which(methods%in%c("rlasso1"))]<-"naive-Lasso-lambda.1se"
            methods[which(methods%in%c("rada"))]<-"naive-AdaLas-lambda.min"
            methods[which(methods%in%c("rada1"))]<-"naive-AdaLas-lambda.1se"
            methods[which(methods%in%c("rX"))]<-"OLS"
            list("m1"=methods,"m2"=methods2)
        }
        cmrr<-cm(colnames(rrid3))
        cmbg<-cm(colnames(bgid3))
        print(colnames(rrid3))
        #print(rrid3)
        if(!is.null(skipid)){
            print("skip")
            df<-data.frame(SqrtSize=rep(sqrt(n_int_list)[-skipid],ncol(rrid3)),
                           Method=rep(cmrr$m1,1,each=nrow(rrid3)),
                           DataSource=rep(cmrr$m2,1,each=nrow(rrid3)),
                           R2=c(matrix(rrid3,nrow=1,byrow = T)))
        }else{
            df<-data.frame(SqrtSize=rep(sqrt(n_int_list),ncol(rrid3)),
                           Method=rep(cmrr$m1,1,each=nrow(rrid3)),
                           DataSource=rep(cmrr$m2,1,each=nrow(rrid3)),
                           R2=c(matrix(rrid3,nrow=1,byrow = T)))
            df2<-data.frame(SqrtSize=rep(sqrt(n_int_list),ncol(bgid3)),
                            Method=rep(cmbg$m1,1,each=nrow(bgid3)),
                            DataSource=rep(cmbg$m2,1,each=nrow(bgid3)),
                            BetaMSE=c(matrix(bgid3,nrow=1,byrow = T)))
        }

        library(ggplot2)
        gg<-ggplot(data = df)+
            geom_point(aes(x=SqrtSize,y=R2,color=Method,shape = DataSource))+
            geom_line(aes(x=SqrtSize,y=R2,color=Method,linetype = DataSource))+
            geom_hline(yintercept = highbase,color = "darkred")+
            annotate("text", x = 20, y = highbase*1.015, label = "Max Achievable R2",size=3)+
            labs(title = paste0("R2:Comp:Incomp=1:",inv_r,";RiskFactor:",p_X))+
            theme_bar()
        ggsave(paste0("plotlocal/",gname,":R2:Comp:Incomp=1:",inv_r,";RiskFactor:",p_X,".png"),gg,width = 7,height = 5)
        gg<-ggplot(data = df2)+
            geom_point(aes(x=SqrtSize,y=BetaMSE,color=Method))+
            geom_line(aes(x=SqrtSize,y=BetaMSE,color=Method,linetype = DataSource))+
            labs(title = paste0("BetaMSE:Comp:Incomp=1:",inv_r,";RiskFactor:",p_X))+
            theme_bar()
        ggsave(paste0("plotlocal/",gname,":BetaMSE:Comp:Incomp=1:",inv_r,";RiskFactor:",p_X,".png"),gg,width = 7,height = 5)
    }
}








