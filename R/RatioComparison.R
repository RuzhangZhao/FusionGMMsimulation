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
gname="Lasso"
gname="AdaptiveLasso"
ratios<-c()
for (id1 in c(1:2)){
    p_X<-p_X_list[id1]
    for(id2 in c(1:2)){
        inv_r<-inv_r_list[id2]
        message(c(id1,id2))
        for(id3 in c(1:length(n_int_list))){
            p_X<-p_X_list[id1]
            inv_r = inv_r_list[id2]
            n_int = n_int_list[id3]
            ratios1<-c()
            res<-readRDS(paste0('../FusionGMMdata/ratio/Ln_n_',n_int,'_ir_',inv_r,'_pX_',p_X,'.rds'))
            for(i in 1:length(res)){
                ratios1<-c(ratios1,res[[i]]$ratios)
            }
            ratios<-rbind(ratios,c(mean(ratios1),p_X,inv_r,n_int))
        }
    }
}
for(i in 1:100){
    print(i)
    a<-lm(y~n_int,data = data.frame(y=ratios[,1],inv_r=ratios[,3]^(1/2),n_int=(ratios[,4])^(1/3)))
    print(summary(a)$r.squared)
}

