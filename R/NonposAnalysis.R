df1<-c()
df2<-c()
df3<-c()
ntest_list<-500*2^(0:5)*30
ntest_list<-c(1:6)*10000
for(i in c(1:6)){
  n<-ntest_list[i]
  #nonpos<-readRDS(paste0("Lnfinal_pos_ir_30_pX_40_n_",n,".rds"))
  nonpos<-readRDS(paste0("pos_ir_30_pX_40_n_",n,".rds"))
  aa<-sapply(c(1:5), function(j){
    bb<-(sapply(1:length(nonpos), function(i){
      nonpos[[i]][[j]]}))
    rowMeans(bb)})
  print(aa)
  for(j in 1:5){
    df1<-rbind(df1,c(rownames(aa)[j],aa[j,1],n))
    df2<-rbind(df2,c(rownames(aa)[j],aa[j,3],n))
    df3<-rbind(df3,c(rownames(aa)[j],aa[j,5],n))
  }
}
cn<-function(df1){
  df1$TLGMM[df1$TLGMM == 'mul_ada']<-'multivariate-2lam'
  df1$TLGMM[df1$TLGMM == 'mul_ada_1lam']<-'multivariate-1lam'
  df1$TLGMM[df1$TLGMM == 'uni_ada']<-'univariate-2lam'
  df1$TLGMM[df1$TLGMM == 'uni_ada_1lam']<-'univariate-1lam'
  df1
}
df1<-data.frame(df1)
colnames(df1)<-c("TLGMM","FDR","n")
df1$n<-round(as.numeric(df1$n)/30)
df1$n<-factor(df1$n,levels = sort(unique(df1$n)))
df1$FDR<-as.numeric(df1$FDR)
df1<-cn(df1)
library(ggplot2)
# Basic barplot
p<-ggplot(data=df1, aes(x=n, y=FDR,fill=TLGMM)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  coord_cartesian(ylim= c(0,0.0505))+
  #labs(title = "Logistic: False Discovery Rate: TLGMM-adaptivelasso")+
  labs(title = "Linear: False Discovery Rate: TLGMM-adaptivelasso")+
  xlab("Internal Sample Size")+ylab("False Discovery Rate")+
  geom_bar(stat="identity", position=position_dodge())+theme_bar()
ggsave(paste0("/users/rzhao1/intgmm/plot/ln_fdr.pdf"),p,height = 4,width = 7)



df2<-data.frame(df2)
colnames(df2)<-c("TLGMM","Power","n")
df2$n<-round(as.numeric(df2$n)/30)
df2$Power<-as.numeric(df2$Power)
df2$n<-factor(df2$n,levels = sort(unique(df2$n)))
df2<-cn(df2)
p<-ggplot(data=df2, aes(x=n, y=Power,fill=TLGMM)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
  coord_cartesian(ylim = c(0,1))+
  #labs(title = "Logistic: Power: TLGMM-adaptivelasso")+
  labs(title = "Linear: Power: TLGMM-adaptivelasso")+
  xlab("Internal Sample Size")+ylab("Power")+
  geom_bar(stat="identity", position=position_dodge())+theme_bar()
ggsave(paste0("/users/rzhao1/intgmm/plot/ln_power.pdf"),p,height = 4,width = 7)



df3<-data.frame(df3)
colnames(df3)<-c("TLGMM","CoverP","n")
df3$n<-round(as.numeric(df3$n)/30)
df3$CoverP<-as.numeric(df3$CoverP)
df3$n<-factor(df3$n,levels = sort(unique(df3$n)))
df3<-cn(df3)
p<-ggplot(data=df3, aes(x=n, y=CoverP,fill=TLGMM)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
  coord_cartesian(ylim = c(0.55,1))+
  #labs(title = "Logistic: Coverage Probability: TLGMM-adaptivelasso")+
  labs(title = "Linear: Coverage Probability: TLGMM-adaptivelasso")+
  xlab("Internal Sample Size")+ylab("Coverage Probability")+
  geom_bar(stat="identity", position=position_dodge())+theme_bar()
ggsave(paste0("/users/rzhao1/intgmm/plot/ln_cover.pdf"),p,height = 4,width = 7)



alpha_list<-c(0.01,0.05,0.1)
p_X<-40
p_A<-150
coef_X<-rep(0,p_X)
coef_A<-rep(0,p_A)
X_nonzero_index<-1:10
A_nonzero_index<-c(11,12,13,17,21,26,31,34,36,50,64,67,69,88,90)-10 #15 nonzero
coef_X[X_nonzero_index]<- 0.3
coef_A[A_nonzero_index]<- 0.25
coefXA<-c(coef_X,coef_A)
#correct_pos<-c(which(coefXA!=0),length(coefXA)+1)
correct_pos<-c(which(coefXA!=0))
method_alpha<-c(0,0,0,0,0)
for(n in 2:6){
  n<-n*10000
  nonpos<-readRDS(paste0("Lg2_pos_ir_30_pX_40_n_",n,".rds"))
  for(j in 1:5){
    aa<-sapply(1:length(nonpos), function(i){
      alpha_cutoff = 0.05
      pval<-nonpos[[i]][[6]][[j]][nonpos[[i]][[7]][[j]]%in%c(1:190)]
      pval1<-p.adjust(pval,method = "BH")
      #pval1<-pval
      pos<-nonpos[[i]][[7]][[j]][nonpos[[i]][[7]][[j]]%in%c(1:190)]
      #print(c(round(min(pval1),3),alpha_cutoff,length(pos[pval1<alpha_cutoff])))
      sum(!pos[pval1<alpha_cutoff]%in%correct_pos)/max(1e-10,length(pos[pval1<alpha_cutoff]))
      #sum(pos[pval1<alpha_cutoff]%in%correct_pos)/length(correct_pos)
      #print(pos[pval1<alpha_cutoff])
    })
    method_alpha[j]<-mean(aa)
  }
  print(method_alpha)
}


for(n in 2:6){
  n<-n*10000
  nonpos<-readRDS(paste0("Lg2_pos_ir_30_pX_40_n_",n,".rds"))
  for(j in 1:5){
    aa<-sapply(1:length(nonpos), function(i){
      alpha_cutoff = 0.05
      pval<-nonpos[[i]][[6]][[j]][nonpos[[i]][[7]][[j]]%in%c(1:190)]
      pval1<-p.adjust(pval,method = "none")
      #pval1<-pval
      pos<-nonpos[[i]][[7]][[j]][nonpos[[i]][[7]][[j]]%in%c(1:190)]
      #print(c(round(min(pval1),3),alpha_cutoff,length(pos[pval1<alpha_cutoff])))
      sum(!pos[pval1<alpha_cutoff]%in%correct_pos)/max(1e-10,length(pos[pval1<alpha_cutoff]))
      #sum(pos[pval1<alpha_cutoff]%in%correct_pos)/length(correct_pos)
    })
    method_alpha[j]<-mean(aa)
  }
  print(method_alpha)
}

namelist<-c("ada",'mul_ada', 'mul_ada_1lam','uni_ada','uni_ada_1lam')
df1<-c()
for(n in c(2:6)){
  n<-n*10000
  nonpos<-readRDS(paste0("Lg2_pos_ir_30_pX_40_n_",n,".rds"))
  for(j in 2:5){
  aa<-aa<-sapply(1:length(nonpos), function(i){
    alpha_cutoff = 0.05
    pval<-nonpos[[i]][[6]][[j]][nonpos[[i]][[7]][[j]]%in%c(1:190)]
    pval1<-p.adjust(pval,method = "BH")
    #pval1<-pval
    pos<-nonpos[[i]][[7]][[j]][nonpos[[i]][[7]][[j]]%in%c(1:190)]
    #print(c(round(min(pval1),3),alpha_cutoff,length(pos[pval1<alpha_cutoff])))
    sum(!pos[pval1<alpha_cutoff]%in%correct_pos)/max(1e-10,length(pos[pval1<alpha_cutoff]))
    #sum(pos[pval1<alpha_cutoff]%in%correct_pos)/length(correct_pos)
  })
    df1<-rbind(df1,c(namelist[j],mean(aa),n))
  }
}

cn<-function(df1){
  df1$TLGMM[df1$TLGMM == 'mul_ada']<-'multivariate-2lam'
  df1$TLGMM[df1$TLGMM == 'mul_ada_1lam']<-'multivariate-1lam'
  df1$TLGMM[df1$TLGMM == 'uni_ada']<-'univariate-2lam'
  df1$TLGMM[df1$TLGMM == 'uni_ada_1lam']<-'univariate-1lam'
  df1
}
df1<-data.frame(df1)
colnames(df1)<-c("TLGMM","FDR","n")
df1$n<-round(as.numeric(df1$n)/30)
df1$n<-factor(df1$n,levels = sort(unique(df1$n)))
df1$FDR<-as.numeric(df1$FDR)
df1<-cn(df1)
library(ggplot2)
# Basic barplot
p<-ggplot(data=df1, aes(x=n, y=FDR,fill=TLGMM)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  coord_cartesian(ylim= c(0.01,0.071))+
  labs(title = "Logistic: False Discovery Rate: TLGMM-adaptivelasso")+
  xlab("Internal Sample Size")+ylab("False Discovery Rate")+
  geom_bar(stat="identity", position=position_dodge())+theme_bar()
ggsave(paste0("/users/rzhao1/intgmm/plot/lg_fdr.pdf"),p,height = 4,width = 7)













df1<-c()
df3<-c()
for(n in c(2:6)){
  n<-n*10000
  nonpos<-readRDS(paste0("Lg3_pos_ir_30_pX_40_n_",n,".rds"))
  aa<-sapply(c(2), function(j){
    bb<-(sapply(1:length(nonpos), function(i){
      nonpos[[i]][[j]]}))
    rowMeans(bb)})
  aa<-c(round(n/30,0),aa)
  names(aa)<-c('internal_size','adalasso','uni_ada_2lam','uni_ada_1lam','mul_ada_2lam','mul_ada_1lam')
  print(aa)
  #for(j in 1:5){
  #  df1<-rbind(df1,c(rownames(aa)[j],aa[j,1],n))
  #  df3<-rbind(df3,c(rownames(aa)[j],aa[j,3],n))
  #}
}
cn<-function(df1){
  df1$TLGMM[df1$TLGMM == 'mul_ada']<-'multivariate-2lam'
  df1$TLGMM[df1$TLGMM == 'mul_ada_1lam']<-'multivariate-1lam'
  df1$TLGMM[df1$TLGMM == 'uni_ada']<-'univariate-2lam'
  df1$TLGMM[df1$TLGMM == 'uni_ada_1lam']<-'univariate-1lam'
  df1
}
df1<-data.frame(df1)
colnames(df1)<-c("TLGMM","FDR","n")
df1$n<-round(as.numeric(df1$n)/30)
df1$n<-factor(df1$n,levels = sort(unique(df1$n)))
df1$FDR<-as.numeric(df1$FDR)
df1<-cn(df1)

library(ggplot2)
# Basic barplot
p<-ggplot(data=df1, aes(x=n, y=FDR,fill=TLGMM)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  coord_cartesian(ylim= c(0.01,0.0725))+
  labs(title = "Logistic: False Discovery Rate: TLGMM-adaptivelasso")+
  xlab("Internal Sample Size")+ylab("False Discovery Rate")+
  geom_bar(stat="identity", position=position_dodge())+theme_bar()
ggsave(paste0("/users/rzhao1/TLGMM/plot/lg_fdr.pdf"),p,height = 4,width = 7)


df3<-data.frame(df3)
colnames(df3)<-c("TLGMM","CoverP","n")
df3$n<-round(as.numeric(df3$n)/30)
df3$CoverP<-as.numeric(df3$CoverP)
df3$n<-factor(df3$n,levels = sort(unique(df3$n)))
df3<-cn(df3)
p<-ggplot(data=df3, aes(x=n, y=CoverP,fill=TLGMM)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
  coord_cartesian(ylim = c(0.75,1))+
  labs(title = "Logistic: Coverage Probability: TLGMM-adaptivelasso")+
  xlab("Internal Sample Size")+ylab("Coverage Probability")+
  geom_bar(stat="identity", position=position_dodge())+theme_bar()
ggsave(paste0("/users/rzhao1/TLGMM/plot/lg_cover.pdf"),p,height = 4,width = 7)


