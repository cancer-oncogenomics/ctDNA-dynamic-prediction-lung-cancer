## Run all of the analyses
library(JMbayes)
library(dplyr)
library("parallel")
library(ggplot2)
require(survival)
require(survminer)
library(caret)
require(R.utils)
library(reshape2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("function/cross_validation.R")
source("function/plot_evaluation.R")
source("function/leaveoneout.R")
source("function/reliability_diagram.R")

#### ---- Run all analysese scripts ---- ####

mydata = read.csv("data/Joint_model_data.csv") %>% mutate(T_stage = as.factor(T_stage),TP53 = as.factor(TP53))
mydata$logAF = log(mydata$meanAF+1e-6)-log(1e-6)

mydata.id = mydata[which(!duplicated(mydata$PatientID)), ]
mydata.id.AD =  subset(mydata.id,PathologicalType == 'AD')


Tstart = 244/30 # landmarking time, using data up to plasma time point 4
predict_times = c(12, 15) # horizon time 

##### 5-fold cross validation for 20 times #####
# set.seed(888)
# splits=createMultiFolds(factor(mydata.id$DFS_status),k=5,times=5)
# splits_AD=createMultiFolds(factor(mydata.id.AD$DFS_status),k=5,times=10)

## For reproduction, we are using cross validation splits used in the paper
splits = readRDS('data/splits.RDS')
splits_AD = readRDS('data/splits.AD.RDS')

### run cross validation for models using all patients

## using multiple cores
cores = parallel::detectCores()
cl <- makeCluster(10)
clusterExport(cl,c("mydata","mydata.id","runCV",'Tstart','predict_times','mydata.id.AD'))
clusterEvalQ(cl, {library(splines);library(nlme);library("JMbayes");
  library(dplyr);library(survival);library("xtable");
  library("lattice")})
res_all <- parLapply(cl, splits, function(x) runCV(x,Tstart,predict_times,mydata,mydata.id))
res_AD <- parLapply(cl, splits_AD, function(x) runCV(x,Tstart,predict_times,mydata,mydata.id.AD))
stopCluster(cl)


# result 
out_all = bind_rows(res_all, .id = "column_label")
out_AD = bind_rows(res_AD, .id = "column_label")

## plot the CV results
                    
# CV_all_patients/JMvsCox_testing.pdf ~ Fig 4b,  CV_all_patients/JMvsCox_training.pdf ~ Supplementary Fig 10,
# CV_all_patients/betweenJMs.pdf ~ Supplementary Fig 15               
plot_evaluation('results/CV_all_patients',out_all) 

# CV_AD_patients/JMvsCox_testing.pdf,JMvsCox_training.pdf ~ Supplementary Fig 12
plot_evaluation('results/CV_AD_patients',out_AD)

## plot the personalized prediction ~ Fig 4c,d; Supplementary Fig 13
lm <- lme( logAF ~ ns(TestDate,2),data=mydata,random = ~  ns(TestDate,2) | PatientID,
           control = lmeControl(opt = "optim",msMaxIter =1000))
fit <- coxph(Surv(DFS,DFS_status)~TP53+T_stage ,data=mydata.id,x = TRUE)
iForm <- list(fixed = ~ 0 + TestDate + ins(TestDate, 2), random = ~ 0 + TestDate + ins(TestDate, 2),
              indFixed = 1:3, indRandom = 1:3)
final_jm = jointModelBayes(lm, fit, timeVar = "TestDate",
                param = "td-extra", extraForm = iForm)

if (!file.exists('personalized')) dir.create('personalized')

for (PID in unique(mydata.id$PatientID)){
  
  ND = mydata[mydata$PatientID==PID,]
  l = nrow(ND)
  Relapse_status = ifelse(ND$DFS_status[1]==1,'Relapsed','Relaspe-free')
  if (l<2) next
  
  png(paste0("personalized/",PID,".png"),width=(3*(l-1)+1)*100,height = 4*100,pointsize = 12,bg='transparent')
  
  survPreds <- vector("list", nrow(ND))
  for (i in 1:nrow(ND)) {
    Tstart = ND[i,"TestDate"]
    survPreds[[i]] <- survfitJM(final_jm,idVar = "PatientID", newdata = ND[1:i, ],
                                survTimes= seq(Tstart,min(max(Tstart+180/30,540/30),570/30),10/30),
    )
  }
  par(mfrow = c(1, l-1),oma = c(0, 2, 2, 2)) 
  for ( i in c(2:l)){
    plot(survPreds[[i]], estimator = "median",include.y = T,
         main=paste0("Follow-up time(months): ",round(survPreds[[i]]$last.time, 1)),
         xlab = "Time (months)",conf.int = TRUE, ylab = "", ylab2 = "" ,cex.lab =1.5,cex.main=1.5
    )
  }
  mtext("log mean VAF", side = 2, line = -1, outer = TRUE,cex=1.5)
  mtext("Recurrence-free Probability", side = 4, line = -1, outer = TRUE,cex=1.5)
  mtext(paste0(PID,", ",Relapse_status), outer = TRUE, cex = 1.5,side = 3,adj = 0)
  dev.off()
}                    
                    
##### leave-one-out cross-validation #####
## joint model
cores = parallel::detectCores()
cl <- makeCluster(10)
clusterExport(cl,c("mydata","mydata.id","runLOO_jm",'Tstart','predict_times'))
clusterEvalQ(cl, {library(splines);library(nlme);library("JMbayes");
  library(dplyr);library(survival);library("xtable");
  library("lattice")})
res_LOO_jm = parLapply(cl, c(1:nrow(mydata.id)), function(x) runLOO_jm(x,Tstart,predict_times,mydata,mydata.id))
stopCluster(cl)
## landmarking cox model
res_LOO_cox = lapply(c(1:nrow(mydata.id)),function(x) runLOO_cox(x,Tstart,predict_times,mydata,mydata.id))

## result data frame
dt_risk_LOO_jm = do.call("rbind",res_LOO_jm)
dt_risk_LOO_cox = do.call("rbind",res_LOO_cox)
dt_risk_LOO = merge(dt_risk_LOO_jm,dt_risk_LOO_cox,by=colnames(mydata.id),all = TRUE)

## plot the reliable diagram ~ Supplementary Fig 11
# at 12 month
reliability_diagram(
  list(subset(dt_risk_LOO,!is.na(prob1)) %>% select(DFS_status,DFS,prob1) ,
       subset(dt_risk_LOO,!is.na(cox1_prob1))%>% select(DFS_status,DFS,cox1_prob1)) ,
  u=12,stat_type = 'C',bins=5,c('Joint model', 'Cox model'),c('blue', 'red'),"12 Months"
)

# at 15 month
reliability_diagram(
  list(subset(dt_risk_LOO,!is.na(prob2)) %>% select(DFS_status,DFS,prob2) ,
       subset(dt_risk_LOO,!is.na(cox1_prob2))%>% select(DFS_status,DFS,cox1_prob2)) ,
  u=15,stat_type = 'C',bins=5,c('Joint model', 'Cox model'),c('blue', 'red'),"15 Months"
)
