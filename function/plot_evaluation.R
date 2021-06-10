library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(rstatix)
######### test ############
# out = out_all
# folder = 'results/CV_all_patients'
############# res ###############
plot_evaluation = function(folder,out){
  
  if (!file.exists(folder)) dir.create(folder)
  

  
  # out 
  out2 = out %>% rowwise() %>% mutate(repeats = unlist(strsplit(column_label, split="[.]"))[2]) %>% as.data.frame()
  get_out = function(rep,model,time,v,tp,k){
    dt = out2 %>% subset(.,repeats==rep&model2==model&variable==v&type==tp)
    n= length(na.omit(dt[,as.character(time)]))
    mean = mean(dt[,as.character(time)],na.rm = TRUE)
    return(dynGet(k))
  }
  new_out = expand.grid(
    'rep' = unique(out2$rep),
    "model2" = unique(out2$model2) , 
    "time" = as.character(predict_times),
    "variable"= c("AUC","PE"),
    "type"= c('training','testing')
  )
  
  new_out$value = mapply(get_out,rep=new_out$rep, model=new_out$model2,time=new_out$time,
                         v=new_out$variable,tp=new_out$type,MoreArgs = list(k="mean"))
  NA_rep = new_out[is.na(new_out$value),c('rep','model2','type','time')]
  for ( i in 1:nrow(NA_rep)){
    new_out[new_out$rep==NA_rep$rep[i]&
              new_out$model2==NA_rep$model2[i]&
              new_out$type==NA_rep$type[i]&
              new_out$time==NA_rep$time[i],'value'] = NA
  }
  
  write.csv(new_out %>% mutate(model2 = recode(model2,"LM_bi" = "cox(Landmark)","LM_P2"="cox(Postsurgical)")),
                               file = paste0(folder,"/outdata.csv"),quote = F)
  
  # res
  get_res = function(model,t,v,tp,k){
    dt = new_out %>% subset(.,model2==model&variable==v&type==tp)
    n= length(na.omit(dt[dt$time==t,]))
    mean = mean(dt[dt$time==t,'value'],na.rm = TRUE)
    return(dynGet(k))
  }
  
  res_table = expand.grid(
    "model2" = unique(new_out$model2), 
    "variable"= c("AUC","PE"),
    "time" = as.character(predict_times),
    "type"= c('training','testing')
  )
  
  res_table$value = mapply(get_res,model=res_table$model2,t=res_table$time,
                           v=res_table$variable,tp=res_table$type,MoreArgs = list(k="mean"))
  res_table = dcast(res_table,...~model2) %>% rename( "cox-Landmark" = LM_bi, 'cox-Postsurgical' = LM_P2) %>% arrange(variable,time,type)
  write.csv(res_table,file = paste0(folder,"/res_table.csv"),quote = F)
  
  
  # label
  rs_label = c("AUC"="AUROC","PE"="Prediction Error")
  rs_label2 = c('12' = "12 months",'15' = "15 months")
  
  ## between JMs using different association structures
  df2 = subset(new_out,model2 %in% c('JM_cumulative','JM_value','JM_value-slope'))  %>%
    rename( evaluation = variable  )
  
  ggplot(df2,aes(x=model2,y=value,fill=type))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.1))+
    scale_fill_brewer(palette ='Pastel1',)+
    scale_color_brewer(palette ='Pastel1')+
    scale_y_continuous(name = "")+
    scale_x_discrete(name="Joint Models",
                     limits=c('JM_value','JM_value-slope','JM_cumulative'),
                     breaks=c('JM_value','JM_value-slope','JM_cumulative'),
                     labels=c('value','value\n+slope','value\n+cumulative'))+
    facet_grid(rows = vars(evaluation) ,cols = vars(time),scales="free",space="free",switch = 'y',
               labeller = labeller(evaluation=rs_label,time = rs_label2))+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      # legend.position = 'none',
      axis.text = element_text(size=12),
      axis.title = element_text(size=15),
      strip.background=element_blank(),
      strip.text=element_text(size = 14),
      strip.placement = 'outside'
    )
  
  ggsave(paste0(folder,"/betweenJMs.pdf"),device='pdf',width=8,height=8)
  
  ## comparison between jm and cox 
  
  df = subset(new_out,model2 %in% c('JM_cumulative','LM_P2','LM_bi')) %>% rename( evaluation = variable  ) 
  
  stat.test <- df  %>%
    group_by(time,evaluation,type) %>%
    wilcox_test(value ~ model2,ref.group = 'JM_cumulative') %>%
    add_significance("p",cutpoints = c( 0,0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>%
    add_xy_position(x='model2',step.increase=0.03)
  
  # df  %>% group_by(time,evaluation,type) %>% wilcox_effsize(value ~ model2,ref.group = 'JM_cumulative')
  
  ggplot(subset(df,type == "testing"),aes(x=model2,y=value,fill=time))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(position = position_jitterdodge(jitter.width=0.3),shape = 21)+
    scale_fill_brewer(palette ='Pastel1')+
    scale_color_brewer(palette ='Pastel1')+
    scale_y_continuous(name = "")+
    scale_x_discrete(name="Models",
                     breaks=c('JM_cumulative','LM_P2','LM_bi'),
                     labels=c('Joint model','Postsurgical\nCox','Landmark\nCox'))+
    # labs(title = paste0('Testing Sets'))+
    facet_grid(rows = vars(evaluation) ,cols = vars(time),scales="free",space="free",switch = 'y',
               labeller = labeller(evaluation=rs_label,time = rs_label2))+
    stat_pvalue_manual(subset(stat.test,type == "testing"), label = "p.signif",hide.ns = F,tip.length = 0)+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      legend.position = 'none',
      axis.text = element_text(size=10),
      axis.title = element_text(size=12),
      strip.background=element_blank(),
      strip.text=element_text(size = 12),
      strip.placement = 'outside'
    )
  
  ggsave(paste0(folder,"/JMvsCox_testing.pdf"),device='pdf',width=6,height=6)
  
  ggplot(subset(df,type == "training"),aes(x=model2,y=value,fill=time))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.3))+
    scale_fill_brewer(palette ='Pastel1')+
    scale_color_brewer(palette ='Pastel1')+
    scale_y_continuous(name = "")+
    scale_x_discrete(name="Models",
                     breaks=c('JM_cumulative','LM_P2','LM_bi'),
                     labels=c('Joint model','Postsurgical\nCox','Landmark\nCox'))+
    # labs(title = paste0('Training Sets'))+
    facet_grid(rows = vars(evaluation) ,cols = vars(time),scales="free",space="free",switch = 'y',
               labeller = labeller(evaluation=rs_label,time = rs_label2))+
    stat_pvalue_manual(subset(stat.test,type == "testing"), label = "p.signif",hide.ns = F,tip.length = 0)+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      legend.position = 'none',
      axis.text = element_text(size=12),
      axis.title = element_text(size=15),
      strip.background=element_blank(),
      strip.text=element_text(size = 14),
      strip.placement = 'new_outside'
    )
  ggsave(paste0(folder,"/JMvsCox_training.pdf"),device='pdf',width=6,height=6)
}

