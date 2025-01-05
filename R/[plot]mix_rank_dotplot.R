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

Method_names_base<-c("random","mv_ct","mv_nc","mv_lognc_scran(1mv)*","mv_PFlogPF",
                     "logmv_ct_seuratv3*","logmv_nc","logmv_lognc(2lmv)","logmv_PFlogPF",
                     "poisson_ct_scran(3pos)*",
                     "disp_ct","disp_nc_seuratv1(4dis)*","mvp_nc_seuratv2*","disp_lognc","disp_PFlogPF",
                     "scanpy_cell_ranger*","mean_max_ct","mean_max_nc(5max)","mean_max_lognc","mean_max_PFlogPF","SCT*")
type0 = "mixplot"
names_list<-c("1mv","2lmv","3pos","4dis","5max")

num_list<-list(c(1,2),c(1,3),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5),c(3,4),c(3,5),c(4,5),
               c(1,2,3),c(1,2,4),c(1,2,5),c(1,3,4),c(1,3,5),c(1,4,5),
               c(2,3,4),c(2,3,5),c(2,4,5),c(3,4,5),
               c(1,2,3,4),c(1,2,3,5),c(1,2,4,5),c(1,3,4,5),c(2,3,4,5),
               c(1,2,3,4,5))
comb_method_names<-c()

for( i in 1:length(num_list)){
    comb_method_names<-c(comb_method_names,paste0(names_list[num_list[[i]]],collapse = ""))
}

###### cite seq

path = "final_base"
path2 = "final_mix"
method_count = 47
#path = "/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/"
data_rank_CITE_all<-matrix(0,nrow = 6,ncol=method_count)
data_rank_multiome_all<-matrix(0,nrow = 6,ncol=method_count)
data_rank_sort_all<-matrix(0,nrow = 6,ncol=method_count)

for( dataset_name in c("bmcite", "cbmc8k","cbmc_pbmc","fetalBM","CD34","seurat_cite","Sucovid"
)){
    print(dataset_name)
    dataset_name0 = ""
    res<-readRDS(paste0(path,dataset_name0,"/final_",dataset_name,"_evaluate_base.rds"))
    res2<-readRDS(paste0(path2,dataset_name0,"/final_",dataset_name,"_evaluate_mix.rds"))
    newinv_knn_ratio<-c(1/res$knn_ratio,1/res2$knn_ratio)
    newinv_nnacc<-c(1/res$`3nn`,1/res2$`3nn`)
    newdist_cor<-c(res$dist_cor,res2$dist_cor)
    
    newinv_var_ratio<-c(1/res$max_var_ratio,1/res2$max_var_ratio)
    newari = c(res$max_ari,res2$max_ari)
    newnmi = c(res$max_nmi,res2$max_nmi)
    
    get_rank<-function(ratio,dec=F) match(c(1:length(ratio)),order(ratio,decreasing = dec))
    rank_table<-rbind(
        get_rank(c(newinv_var_ratio),dec = T),
        get_rank(c(newinv_knn_ratio),dec = T),
        get_rank(c(newinv_nnacc),dec = T),
        get_rank(c(newdist_cor),dec = T),
        get_rank(c(newari),dec = T),
        get_rank(c(newnmi),dec = T)
    )
    data_rank_CITE_all<-data_rank_CITE_all+rank_table
}

for( dataset_name in c("pbmc3k","human_brain_3k","mouse_brain_fresh_5k","pbmc10k","lymphoma_14k"
)){
    dataset_name0 = ""
    print(dataset_name)
    res<-readRDS(paste0(path,dataset_name0,"/final_",dataset_name,"_evaluate_base.rds"))
    res2<-readRDS(paste0(path2,dataset_name0,"/final_",dataset_name,"_evaluate_mix.rds"))
    newm_inv_knn_ratio<-c(1/res$knn_ratio,1/res2$knn_ratio)
    newm_inv_var_ratio<-c(1/res$var_ratio,1/res2$var_ratio)
    newm_inv_nnacc<-c(1/res$`3nn`,1/res2$`3nn`)
    newm_dist_cor<-c(res$dist_cor,res2$dist_cor)
    newm_ari = c(res$ari,res2$ari)
    newm_nmi = c(res$nmi,res2$nmi)
    
    newm_inv_var_ratio<-c(1/res$max_var_ratio,1/res2$max_var_ratio)
    newm_ari = c(res$max_ari,res2$max_ari)
    print("ari")
    print(newm_ari)
    newm_nmi = c(res$max_nmi,res2$max_nmi)
    print("nmi")
    print(newm_nmi)
    
    rank_table<-rbind(
        get_rank(c(newm_inv_var_ratio),dec = T),
        get_rank(c(newm_inv_knn_ratio),dec = T),
        get_rank(c(newm_inv_nnacc),dec = T),
        get_rank(c(newm_dist_cor),dec = T),
        get_rank(c(newm_ari),dec = T),
        get_rank(c(newm_nmi),dec = T)#,
    )
    
    data_rank_multiome_all<-data_rank_multiome_all+rank_table
}

for( dataset_name in 
     c("duo8","duo4un","duo4","human","mouse","GBM_sd","zheng")){
    dataset_name0 = ""
    
    res<-readRDS(paste0(path,dataset_name0,"/final_",dataset_name,"_evaluate_base.rds")) 
    
    res2<-readRDS(paste0(path2,dataset_name0,"/final_",dataset_name,"_evaluate_mix.rds"))
    
    newsinv_var_ratio<-1/c(res$var_ratio,res2$var_ratio)
    newsnnacc<-c(res$`3nn`,res2$`3nn`)
    newsdist_cor<-1/c(res$dist_ratio,res2$dist_ratio)
    newlisi<-1/c(res$lisi,res2$lisi)
    newasw<-c(res$asw_score,res2$asw_score)
    newsari<-c(res$ari,res2$ari)
    newnmi<-c(res$nmi,res2$nmi)
    newf1<-c(res$f1,res2$f1)
    
    newsari<-c(res$max_ari,res2$max_ari)
    newnmi<-c(res$max_nmi,res2$max_nmi)
    rank_table<-rbind(
        get_rank(c(newsinv_var_ratio),dec=T),
        get_rank(c(newsari),dec = T),
        get_rank(c(newsnnacc),dec=T),
        get_rank(newlisi,dec=T),
        get_rank(newasw,dec=T),
        get_rank(newnmi,dec=T)#,
    )
    
    data_rank_sort_all<-data_rank_sort_all+rank_table
}


library(ggplot2)
Method_names<-c(Method_names_base,comb_method_names)
MethodType1<-c("random",
               rep("mean-var(mv)",4),
               rep("logmean-logvar(logmv)",4),
               "poisson_scran",
               rep("dispersion(disp)",6),
               rep("mean_max",4),
               "SCT")

MethodType2<-c(rep("2Mix",10),rep("3Mix",10),
               rep("4Mix",5),rep("5Mix",1))
MethodType<-c(MethodType1,MethodType2)

ave_rank<-colMeans(rbind(data_rank_CITE_all,data_rank_multiome_all,data_rank_sort_all))

ave_rank<-max(ave_rank)/ave_rank
csdn_bar<-data.frame(Method=Method_names,Rank=ave_rank)
csdn_bar$Method<-factor(csdn_bar$Method,levels = Method_names)
csdn_bar$MethodType<-factor(MethodType,levels = unique(MethodType))
q1<-ggplot(data=csdn_bar, mapping=aes(x = Method, y =Rank, fill=MethodType))+
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6)+
    ylab('') +
    theme(axis.ticks = element_blank())+
    theme_bar()+
    theme(axis.text.x = element_text(angle = 55, hjust = 1))+
    theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+scale_fill_manual(values =c("#808080","#00AFBB","#1B9E77","#FC074E","#4682B4","#E7B800","#BEAED4","#FDC086","#D95F02","#7570B3","#66A61E","#FC4E07"))
ggsave(paste0("plot_finalmix/oraverank_",type0,".pdf"),q1,width = 10,height = 5)


for(i in 1:nrow(data_rank_CITE_all)){
    data_rank_CITE_all[i,]<-max(data_rank_CITE_all[i,])/data_rank_CITE_all[i,]
}
for(i in 1:nrow(data_rank_multiome_all)){
    data_rank_multiome_all[i,]<-max(data_rank_multiome_all[i,])/data_rank_multiome_all[i,] 
}
for(i in 1:nrow(data_rank_sort_all)){
    data_rank_sort_all[i,]<-max(data_rank_sort_all[i,])/data_rank_sort_all[i,]
}



### Allcomb 

rank_table<-rbind(data_rank_sort_all,data_rank_CITE_all,data_rank_multiome_all)
part_index=1:nrow(rank_table)
rank_table<-rank_table[part_index,]

crit = c("1/var_ratio[sort]","ari[sort]",
         "knn[sort]",
         "1/lisi[sort]","asw[sort]","nmi[sort]",
         "1/var_ratio[cite]","1/knn_dist_ratio[cite]",
         "1/knn_mse[cite]","dist_cor[cite]",
         "ari[cite]","nmi[cite]",
         "1/var_ratio[mult]","1/knn_dist_ratio[mult]",
         "1/knn_mse[mult]","dist_cor[mult]",
         "ari[mult]","nmi[mult]")
crit<-crit[part_index]
library(ggplot2)
Method_names<-c(Method_names_base,comb_method_names)
colnames(rank_table) = Method_names
rownames(rank_table) = crit 

rank_df<-c()
for(i in 1:ncol(rank_table)){
    for(j in 1:nrow(rank_table)){
        rank_df<-rbind(rank_df,c(colnames(rank_table)[i],rownames(rank_table)[j],rank_table[j,i],(rank_table[j,i])/(max(rank_table)),-1+exp((rank_table[j,i])/(max(rank_table[j,]))),log(1+(rank_table[j,i])/(max(rank_table[j,])))) )
    }
}

colnames(rank_df)<-c("HVG_Method","Evaluation","Rank","Rank_perc","Rank_color","Rank_color2")

rank_df<-data.frame(rank_df)
rank_df$Rank<-as.numeric(rank_df$Rank)
rank_df$Rank_perc<-as.numeric(rank_df$Rank_perc)
rank_df$Rank_color<-as.numeric(rank_df$Rank_color)
rank_df$Rank_color2<-as.numeric(rank_df$Rank_color2)
rank_df$HVG_Method<-factor(rank_df$HVG_Method,levels = Method_names)
rank_df$Evaluation<-factor(rank_df$Evaluation,levels = crit)
q1<-ggplot(rank_df) 
typelist<-c("mv_lognc_scran(1mv)*","logmv_lognc(2lmv)","poisson_ct_scran(3pos)*","disp_nc_seuratv1(4dis)*","mean_max_nc(5max)","1mv3pos","1mv2lmv3pos","1mv3pos4dis")
for (type in typelist){
    q1<-q1+geom_hline(yintercept =type,linetype = "twodash",color = "red", size = 0.6,alpha=0.4)
}

q1<-q1+
    geom_point(aes(y=HVG_Method, x = Evaluation, color = Rank_color, size = Rank_perc)) + 
    #cowplot::theme_cowplot() + 
    scale_color_viridis_c(name = "Ave_Rank",option = "D")+
    ylab('') +
    theme(axis.ticks = element_blank())+
    theme_bar()+
    scale_y_discrete(expand = c(0.02, 0.02))+
    theme(axis.text.x = element_text(angle = 35, hjust = 1,size=11))+
    theme(axis.text.y = element_text(hjust = 0.5))

ggsave(paste0("plot_finalmix/allcomb_",type0,".pdf"),q1,width = 10,height = 11)






