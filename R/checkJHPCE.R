p_X_list<-c(10,40)
inv_r_list<-c(10,30)
n_int_list<-floor((seq(sqrt(100),sqrt(3000),(sqrt(3000) - sqrt(100))/9))^2)
for(id1 in c(2)){
    for(id2 in c(1:2)){
        for(id3 in c(3:10)){
            p_X<-p_X_list[id1]
            inv_r<-inv_r_list[id2]
            n_int<-n_int_list[id3]

            a<-readRDS(paste0('/users/rzhao1/intgmm/result/Lg/Lg_n_',n_int,'_ir_',inv_r,'_pX_',p_X,'.rds'))
            if(length(a)==50){
                print(c(id1,id2,id3))
                print(paste0("p_X=",p_X,";inv_r=",inv_r,":n_int=",n_int))
                print(length(a))
            }
            if(length(a)!=50){
                print("WarningWarning!!!!!")
                print(c(id1,id2,id3))
                print(paste0("p_X=",p_X,";inv_r=",inv_r,":n_int=",n_int))
                print(length(a))
            }
            }
        }}
