getsurvplot <- function(feature,feature_name,df,labs=c("ctDNA Negative","ctDNA Positive")){
  s_fit = as.formula(paste0('Surv(DFS,DFS_status)~', feature))
  fit = surv_fit(s_fit,data=df)
  p =  ggsurvplot(fit, risk.table = TRUE,legend.labs=labs,risk.table.title = "No. at risk",
                  palette = c("#003B63","#9E333E"),legend = 'none',
                  title=feature_name,xlab="Days from surgery",ylab="Recurrence-free Survival",
                  conf.int= F,break.time.by = 90,axes.offset=F,fontsize = 5,
                  tables.theme = theme_cleantable(base_family= 'Arial'),
                  risk.table.height = 0.2,censor.size=2, size = 0.3,
                  ggtheme = theme(
                    axis.text = element_text(size=12),
                    panel.background = element_rect(fill = "white"),
                    plot.title = element_text(hjust = 0.5),
                    plot.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.title = element_text(size=14,face = 'bold', family="sans"),
                    legend.text = element_text(size=14, family="sans"),
                    legend.title = element_text(size=14, family="sans"),
                    title = element_text(size=14,face = 'bold', family="sans"),
                    # panel.grid.major.y = element_line(colour = "grey"),
                    #panel.grid.minor.y = element_line(colour = "darkgrey")
                  ))
  res<-coxph(s_fit,data=df)
  pval = summary(res)$logtest[3]
  if (pval<0.0001){
    pcut = "P < .0001"
  }else if (pval<0.001){
    pcut = "P < .001"
  }else if(pval<0.01){
    pcut = "P < .01"
  }else if(pval<0.05){
    pcut= "P < .05"
  }else{
    pcut= paste0("P = ",round(pval,2))
  }
  p$plot = p$plot+
    ggplot2::annotate("text",x = 30, y = 0.10,size=5,hjust = 0,
                      label = paste0("Log-rank ",pcut,
                                     "\n",
                                     "HR: ",round(summary(res)$conf.int[1],2)," (",
                                     "95%CI: ",round(summary(res)$conf.int[3],2),"-",
                                     round(summary(res)$conf.int[4],2),")"
                      )
    )
  p$table <- p$table + theme(panel.grid.major.y = element_blank(),axis.text.y = element_text(hjust=0,size=14, family="sans"),
                             plot.title = element_text(hjust = 0,size=14,face = 'bold', family="sans"),
                             text = element_text(size=12, family="sans")
  )
  return(p)
}

getATplot = function(title,AT_df){
  res_pos = summary(coxph(Surv(DFS,DFS_status)~AdjuvantTherapy_status,
                          subset(AT_df, P2_stat=="Pos"&PathologicStage %in% c("IIA","IIB","IIIA","IIIB"))))
  pval_pos = res_pos$sctest[3]
  
  res_neg = summary(coxph(Surv(DFS,DFS_status)~AdjuvantTherapy_status,
                          subset(AT_df, P2_stat=="Neg"&PathologicStage %in% c("IIA","IIB","IIIA","IIIB"))))
  pval_neg = res_neg$sctest[3]
  
  res_ACT = summary(coxph(Surv(DFS,DFS_status)~P2_stat,
                          subset(AT_df, AdjuvantTherapy_status==1&PathologicStage %in% c("IIA","IIB","IIIA","IIIB"))))
  pval_ACT = res_ACT$sctest[3]
  
  res_nonACT = summary(coxph(Surv(DFS,DFS_status)~P2_stat,
                             subset(AT_df, AdjuvantTherapy_status==0&PathologicStage %in% c("IIA","IIB","IIIA","IIIB"))))
  pval_nonACT = res_nonACT$sctest[3]
  
  transp = function(pval){
    if (pval<0.0001){
      pcut = "p < .0001"
    }else if (pval<0.001){
      pcut = "p < .001"
    }else if(pval<0.01){
      pcut = "p < .01"
    }else if(pval<0.05){
      pcut= "p < .05"
    }else{
      pcut= paste0("p = ",round(pval,2))
    }
    return(pcut)
  }
  p = ggsurvplot(survfit(Surv(DFS,DFS_status)~Postsurgical_detection+AdjuvantTherapy_status,
                         subset(AT_df, PathologicStage %in% c("IIA","IIB","IIIA","IIIB"))), 
                 risk.table = TRUE,legend = 'none',risk.table.title = "No. at risk",
                 legend.labs=c("no  neg","yes  neg",
                               "no  pos","yes  pos"),
                 palette = c("cornflowerblue","green","brown4","darkorange2"),
                 xlab="Days from surgery",ylab="Recurrence-free Survival",
                 conf.int= F,break.time.by = 90,axes.offset=F,fontsize = 4,censor.size=2, size = 0.5,
                 tables.theme = theme_cleantable(base_family= 'Courier'),
                 risk.table.height = 0.3,title =title,
                 ggtheme = theme(
                   axis.text = element_text(size=12),
                   panel.background = element_rect(fill = "white"),
                   plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.title = element_text(size=14,face = 'bold'),
                   legend.text = element_text(size=14),
                   legend.title = element_text(size=14),
                   title = element_text(size=14,face = 'bold'),
                   # panel.grid.major.y = element_line(colour = "grey"),
                   #panel.grid.minor.y = element_line(colour = "darkgrey")
                 ))
  
  
  p$plot = p$plot+
    ggplot2::annotate("text",x = 10, y = 0.10,size=5,hjust = 0,vjust=0,
                      label = paste0(
                        "ACT patients, ctDNA pos vs neg: ",transp(pval_ACT) ,"; ",
                        "non-ACT patients, ctDNA pos vs neg: ",transp(pval_nonACT) ,"\n",
                        "ctDNA negative, ACT vs non-ACT: ",transp(pval_neg) ,"; ",
                        "ctDNA positive, ACT vs non-ACT: ",transp(pval_pos)
                        
                      )
    )
  p$table <- p$table + theme(panel.grid.major.y = element_blank(),axis.text.y = element_text(hjust=0,size=14),
                             plot.title = element_text(hjust = 0,size=10,face = 'bold'),
                             text = element_text(size=10)
  )
  p
}
